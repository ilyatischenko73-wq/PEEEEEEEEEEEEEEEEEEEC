"""Новый модуль: H.264/NVDEC -> NVMM -> RF-DETR/CUDA -> GPU OSD -> H.264/NVENC.

Требования: Linux x86_64, NVIDIA dGPU с NVDEC/NVENC, DeepStream + pyds,
PyGObject, CuPy и CUDA-сборка PyTorch. Это не backend OpenCV CAP_GSTREAMER.
Кадр живёт в GstBuffer до завершения pad-probe; DLPack не выносится в очередь.
"""

import ctypes
import logging
import math
import os
from pathlib import Path
import platform
from queue import Full
import re
import signal
import time
from urllib.parse import unquote, urlparse

import torch
from config import RECORD, TEXT_SCALE

logger = logging.getLogger(__name__)


class GPUH264Stream:
    def __init__(self, source, output_url, model, gpu_id, shared_queue, pause_event):
        if platform.system() != "Linux" or platform.machine() != "x86_64":
            raise RuntimeError("Этот NVMM/CuPy мост рассчитан на Linux x86_64 + NVIDIA dGPU.")
        # Новый mux не требует заранее знать размер кадра и не масштабирует видео.
        # Только окружение дочернего worker, до загрузки GStreamer-плагинов.
        os.environ.setdefault("USE_NEW_NVSTREAMMUX", "yes")
        try:
            import gi
            gi.require_version("Gst", "1.0")
            from gi.repository import Gst, GLib
            import pyds
            import cupy as cp
        except (ImportError, ValueError) as exc:
            raise RuntimeError("Нужны DeepStream, соответствующий ему pyds, PyGObject и CuPy; см. README_RU.md.") from exc
        if not hasattr(pyds, "get_nvds_buf_surface_gpu"):
            raise RuntimeError("В установленном pyds нет get_nvds_buf_surface_gpu; нужен совместимый DeepStream Python binding.")

        self.Gst, self.GLib, self.pyds, self.cp = Gst, GLib, pyds, cp
        Gst.init(None)
        self.source, self.output_url = str(source), str(output_url)
        self.model, self.gpu_id = model, int(gpu_id)
        self.shared_queue, self.pause_event = shared_queue, pause_event
        self.is_live = urlparse(self.source).scheme.lower() in ("rtsp", "rtsps")
        if urlparse(self.output_url).scheme.lower() not in ("rtsp", "rtsps"):
            raise ValueError("STREAM_OUT_PATH должен быть RTSP URL сервера с поддержкой публикации.")

        torch.cuda.set_device(self.gpu_id)
        cp.cuda.Device(self.gpu_id).use()
        # Проверяем и конфигурацию RF-DETR, и фактическое устройство весов.
        expected = torch.device(f"cuda:{self.gpu_id}")
        inner = self.model.model.model
        configured = torch.device(inner.device)
        if configured.type == "cuda" and configured.index is None:
            configured = torch.device(f"cuda:{torch.cuda.current_device()}")
        if configured != expected or next(inner.model.parameters()).device != expected:
            raise RuntimeError("RF-DETR и видеоплагины должны использовать один и тот же cuda:N.")

        self.latency_ms = int(os.getenv("GST_RTSP_LATENCY_MS", "100"))
        self.bitrate = int(os.getenv("GST_H264_BITRATE", "3000000"))
        self.timeout = float(os.getenv("GST_STALL_TIMEOUT", "30"))
        if self.latency_ms < 0 or self.bitrate <= 0 or not math.isfinite(self.timeout) or self.timeout <= 0:
            raise ValueError("Некорректные GST_RTSP_LATENCY_MS / GST_H264_BITRATE / GST_STALL_TIMEOUT.")
        self.record_enabled = str(RECORD or "FALSE").strip().lower() == "true"
        self.record_count = 0
        if self.record_enabled:
            # Сохраняем исходную опцию JPEG. Только эта явно включённая ветка
            # скачивает кадр на CPU; обратно в видеотракт он не загружается.
            name = re.sub(r'_+', '_', re.sub(r'[/:*?"<>|\\]', '_', self.source)).strip('_')
            self.record_dir = Path("records") / name
            self.record_dir.mkdir(parents=True, exist_ok=True)
            logger.warning("RECORD=True: JPEG включает GPU->CPU и CPU-кодирование JPEG. Для минимальной нагрузки задайте RECORD=False.")

        self.pipeline = None
        self.loop = GLib.MainLoop()
        self.error = None
        self.restart = False
        self.last_frame = time.monotonic()
        self.last_encoded = self.last_frame
        self.frames = 0
        self.encoded = 0
        self._capsule_pointer = ctypes.pythonapi.PyCapsule_GetPointer
        self._capsule_pointer.restype = ctypes.c_void_p
        self._capsule_pointer.argtypes = [ctypes.py_object, ctypes.c_char_p]

    def _fail(self, error):
        # Исключение Python в pad-probe само по себе не завершает GstPipeline.
        # Передаём ошибку в run(), чтобы существующий NATS monitor перезапустил worker.
        if self.error is None:
            self.error = error
            self.GLib.idle_add(self.loop.quit)

    def _link_source(self, element, pad):
        Gst = self.Gst
        try:
            caps = pad.get_current_caps() or pad.query_caps(None)
            if not caps or caps.is_empty():
                return
            info = caps.get_structure(0)
            kind = info.get_name()
            if self.is_live:
                if kind != "application/x-rtp" or info.get_string("media") != "video":
                    return
                if (info.get_string("encoding-name") or "").upper() != "H264":
                    raise RuntimeError("Входной RTSP должен отдавать H.264. Выберите H.264 в настройках камеры.")
                target = self.pipeline.get_by_name("depay")
            else:
                if not kind.startswith("video/"):
                    return
                if kind != "video/x-h264":
                    raise RuntimeError("Видеодорожка MP4/MOV должна быть H.264.")
                target = self.pipeline.get_by_name("input_parser")
            sink = target.get_static_pad("sink")
            if not sink.is_linked() and pad.link(sink) != Gst.PadLinkReturn.OK:
                raise RuntimeError("Не удалось подключить H.264 дорожку к GStreamer pipeline.")
        except Exception as exc:
            self._fail(exc)

    def _geometry(self, pad, info):
        Gst = self.Gst
        try:
            event = info.get_event()
            if event.type == Gst.EventType.CAPS:
                caps = event.parse_caps()
                structure = caps.get_structure(0)
                if not caps.get_features(0).contains("memory:NVMM"):
                    raise RuntimeError("Декодер вернул CPU memory вместо NVMM.")
                width, height = structure.get_value("width"), structure.get_value("height")
                if width <= 0 or height <= 0 or width % 2 or height % 2:
                    raise RuntimeError("Для H.264/NV12 нужны положительные чётные ширина и высота.")
                logger.info("NVDEC: %sx%s, NVMM, cuda:%s", width, height, self.gpu_id)
        except Exception as exc:
            self._fail(exc)
            return Gst.PadProbeReturn.DROP
        return Gst.PadProbeReturn.OK

    def _build(self):
        Gst = self.Gst
        source_part = (
            "rtspsrc name=source rtph264depay name=depay ! "
            if self.is_live else "filesrc name=source ! qtdemux name=demux "
        )
        # URL задаётся через set_property, а не подставляется в parse_launch.
        # Drop только уже декодированных кадров: ссылки H.264 GOP не повреждаются.
        self.pipeline = Gst.parse_launch(
            source_part +
            "h264parse name=input_parser ! video/x-h264,stream-format=byte-stream,alignment=au ! "
            "nvv4l2decoder name=decoder ! identity name=pace ! "
            "queue max-size-buffers=1 max-size-bytes=0 max-size-time=0 leaky=downstream ! mux.sink_0 "
            "nvstreammux name=mux batch-size=1 batched-push-timeout=33000 ! "
            "nvvideoconvert name=rgba ! video/x-raw(memory:NVMM),format=RGBA ! "
            "identity name=infer ! nvdsosd name=osd ! "
            "nvvideoconvert name=nv12 ! video/x-raw(memory:NVMM),format=NV12 ! "
            "nvv4l2h264enc name=encoder ! video/x-h264,profile=main ! "
            "h264parse name=output_parser config-interval=-1 ! rtspclientsink name=output"
        )
        for name in ("decoder", "rgba", "osd", "nv12", "encoder"):
            self.pipeline.get_by_name(name).set_property("gpu-id", self.gpu_id)
        # В decoder 0 означает DEVICE, а в nvvideoconvert/mux DEVICE = 2.
        self.pipeline.get_by_name("decoder").set_property("cudadec-memtype", 0)
        for name in ("rgba", "nv12"):
            converter = self.pipeline.get_by_name(name)
            converter.set_property("nvbuf-memory-type", 2)
            converter.set_property("compute-hw", 1)
        mux = self.pipeline.get_by_name("mux")
        if mux.find_property("width") is not None:
            raise RuntimeError("Загружен старый nvstreammux. Задайте USE_NEW_NVSTREAMMUX=yes до запуска Python и проверьте DeepStream.")
        for name, value in (("gpu-id", self.gpu_id), ("nvbuf-memory-type", 2), ("live-source", True)):
            if mux.find_property(name) is not None:
                mux.set_property(name, value)
        # Сохраняем PTS, файлы подаём с их временными метками, а не на максимальной скорости.
        self.pipeline.get_by_name("pace").set_property("sync", not self.is_live)
        osd = self.pipeline.get_by_name("osd")
        osd.set_property("process-mode", 1)
        osd.set_property("display-text", True)
        osd.set_property("display-bbox", True)
        encoder = self.pipeline.get_by_name("encoder")
        for name, value in (("bitrate", self.bitrate), ("control-rate", 1), ("iframeinterval", 30),
                            ("idrinterval", 30), ("preset-id", 1), ("tuning-info-id", 3)):
            encoder.set_property(name, value)

        output = self.pipeline.get_by_name("output")
        output.set_property("location", self.output_url)
        output.set_property("protocols", 4)  # GstRtsp.RTSPLowerTrans.TCP
        output.set_property("latency", self.latency_ms)
        output.set_property("tcp-timeout", 5_000_000)  # микросекунды
        source = self.pipeline.get_by_name("source")
        if self.is_live:
            source.set_property("location", self.source)
            source.set_property("protocols", 4)
            source.set_property("latency", self.latency_ms)
            source.set_property("drop-on-latency", True)
            source.set_property("tcp-timeout", 5_000_000)
            source.connect("pad-added", self._link_source)
        else:
            parsed = urlparse(self.source)
            if parsed.scheme and parsed.scheme != "file":
                raise ValueError("Поддерживаются RTSP и локальный MP4/MOV с H.264.")
            path = Path(unquote(parsed.path) if parsed.scheme == "file" else self.source).resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            source.set_property("location", str(path))
            self.pipeline.get_by_name("demux").connect("pad-added", self._link_source)
        self.pipeline.get_by_name("decoder").get_static_pad("src").add_probe(Gst.PadProbeType.EVENT_DOWNSTREAM, self._geometry)
        self.pipeline.get_by_name("infer").get_static_pad("src").add_probe(Gst.PadProbeType.BUFFER, self._process_frame)
        self.pipeline.get_by_name("output_parser").get_static_pad("src").add_probe(Gst.PadProbeType.BUFFER, self._encoded_frame)

    def _display(self, batch, meta, detections, labels, width, height):
        # На CPU передаются только координаты, классы и строки, но не пиксели.
        # nvdsosd рисует по этим метаданным на исходной NVMM-поверхности.
        pyds = self.pyds
        for index, (box, label) in enumerate(zip(detections.xyxy, labels)):
            coords = [float(value) for value in box]
            if not all(math.isfinite(value) for value in coords):
                continue
            x1, y1, x2, y2 = coords
            x1, x2 = max(0.0, min(width, x1)), max(0.0, min(width, x2))
            y1, y2 = max(0.0, min(height, y1)), max(0.0, min(height, y2))
            if x2 <= x1 or y2 <= y1:
                continue
            obj = pyds.nvds_acquire_obj_meta_from_pool(batch)
            if obj is None:
                raise RuntimeError("Не удалось выделить NvDsObjectMeta.")
            obj.object_id = 0xFFFFFFFFFFFFFFFF  # UNTRACKED_OBJECT_ID
            obj.unique_component_id = 1
            obj.class_id = int(detections.class_id[index])
            obj.confidence = float(detections.confidence[index])
            rect = obj.rect_params
            rect.left, rect.top, rect.width, rect.height = x1, y1, x2 - x1, y2 - y1
            rect.border_width = 2
            rect.has_bg_color = 0
            rect.border_color.set(0.0, 1.0, 0.0, 1.0)
            text = obj.text_params
            text.display_text = str(label)
            text.x_offset, text.y_offset = int(x1), max(0, int(y1) - 24)
            text.font_params.font_name = "Serif"
            text.font_params.font_size = max(8, int(32 * float(TEXT_SCALE)))
            text.font_params.font_color.set(1.0, 1.0, 1.0, 1.0)
            text.set_bg_clr = 1
            text.text_bg_clr.set(0.0, 0.0, 0.0, 0.7)
            pyds.nvds_add_obj_meta_to_frame(meta, obj, None)

    def _process_frame(self, pad, info):
        Gst, cp, pyds = self.Gst, self.cp, self.pyds
        if self.error is not None:
            return Gst.PadProbeReturn.DROP
        try:
            self.last_frame = time.monotonic()
            self.frames += 1
            record_this = self.record_enabled and self.record_count % 5 == 0
            self.record_count += 1
            paused = self.pause_event.is_set()
            if paused and not record_this:
                return Gst.PadProbeReturn.OK
            buffer = info.get_buffer()
            batch = pyds.gst_buffer_get_nvds_batch_meta(hash(buffer))
            if batch is None or batch.frame_meta_list is None:
                raise RuntimeError("nvstreammux не добавил NvDsBatchMeta.")
            meta = pyds.NvDsFrameMeta.cast(batch.frame_meta_list.data)
            # batch-size=1: один worker обрабатывает одну камеру.
            with torch.cuda.device(self.gpu_id), cp.cuda.Device(self.gpu_id), torch.inference_mode():
                dtype, shape, strides, capsule, size = pyds.get_nvds_buf_surface_gpu(hash(buffer), meta.batch_id)
                pointer = self._capsule_pointer(capsule, None)
                if not pointer or len(shape) != 3 or shape[2] != 4:
                    raise RuntimeError("Ожидался CUDA-буфер RGBA (H,W,4).")
                # owner=buffer удерживает исходный GstBuffer, strides учитывают pitch.
                memory = cp.cuda.UnownedMemory(pointer, size, buffer, device_id=self.gpu_id)
                rgba_array = cp.ndarray(shape, dtype=dtype, memptr=cp.cuda.MemoryPointer(memory, 0), strides=strides)
                # DeepStream и PyTorch могут использовать разные CUDA streams.
                # Синхронизация не копирует кадр на CPU; завершение до возврата
                # из probe обязательно, иначе decoder сможет переиспользовать память.
                torch.cuda.synchronize(self.gpu_id)
                try:
                    rgba = torch.utils.dlpack.from_dlpack(rgba_array)
                    if rgba.device != torch.device(f"cuda:{self.gpu_id}") or rgba.dtype != torch.uint8:
                        raise RuntimeError("DLPack вернул неправильное CUDA-устройство или dtype.")
                    if record_this:
                        import cv2
                        bgr = rgba[:, :, :3].flip(-1).contiguous().cpu().numpy()
                        filename = self.record_dir / f"{(self.record_count - 1) % 2000}.jpg"
                        if not cv2.imwrite(str(filename), bgr):
                            logger.warning("Не удалось записать JPEG: %s", filename)
                    if not paused:
                        # RGBA -> RGB CHW, float32 [0,1]: только GPU-операции.
                        rgb = rgba[:, :, :3].permute(2, 0, 1).to(dtype=torch.float32).contiguous().div_(255.0)
                        _, detections, labels = self.model.predictImage(rgb, annotate=False, return_labels=True)
                        self._display(batch, meta, detections, labels, shape[1], shape[0])
                        try:
                            self.shared_queue.put_nowait(detections)
                        except Full:
                            pass  # Если очередь ограничена вызывающим кодом, не задерживаем видео.
                finally:
                    torch.cuda.synchronize(self.gpu_id)
        except Exception as exc:
            logger.exception("Ошибка GPU-видеотракта")
            self._fail(exc)
            return Gst.PadProbeReturn.DROP
        return Gst.PadProbeReturn.OK

    def _encoded_frame(self, pad, info):
        self.last_encoded = time.monotonic()
        self.encoded += 1
        if self.encoded == 1:
            logger.info("NVENC выдаёт H.264; GPU-видеотракт запущен.")
        return self.Gst.PadProbeReturn.OK

    def _bus_message(self, bus, message):
        Gst = self.Gst
        if message.type == Gst.MessageType.ERROR:
            error, debug = message.parse_error()
            logger.error("GStreamer %s: %s (%s)", message.src.get_name(), error.message, debug)
            self._fail(RuntimeError(f"GStreamer {message.src.get_name()}: {error.message}"))
        elif message.type == Gst.MessageType.EOS:
            if not self.encoded:
                self._fail(RuntimeError("Источник завершился без выходных H.264 кадров; проверьте видеодорожку и кодеки."))
                return
            # Локальные тестовые MP4, как и раньше, повторяются. Полный restart
            # сбрасывает PTS/сегменты и RTSP-сессию без отправки обратных timestamps.
            self.restart = not self.is_live
            if self.is_live:
                self._fail(RuntimeError("RTSP источник завершился (EOS)."))
            else:
                self.loop.quit()

    def _watchdog(self):
        now = time.monotonic()
        # Даём первому CUDA inference/RTSP handshake больше времени.
        limit = max(60.0, self.timeout) if self.encoded == 0 else self.timeout
        if now - min(self.last_frame, self.last_encoded) > limit:
            self._fail(RuntimeError("Нет новых кадров или H.264 на выходе: проверьте источник, NVDEC/NVENC и RTSP сервер."))
        return True

    def run(self):
        Gst, GLib = self.Gst, self.GLib
        self.stopping = False
        previous_sigterm = signal.getsignal(signal.SIGTERM)
        def stop(signum, frame):
            self.stopping = True
            self.restart = False
            # idle сработает и если SIGTERM пришёл непосредственно перед loop.run().
            GLib.idle_add(self.loop.quit)
        signal.signal(signal.SIGTERM, stop)
        try:
            while True:
                self.restart = False
                self.last_frame = self.last_encoded = time.monotonic()
                self.frames = self.encoded = 0
                bus, handler, timer = None, None, None
                try:
                    self._build()
                    bus = self.pipeline.get_bus()
                    bus.add_signal_watch()
                    handler = bus.connect("message", self._bus_message)
                    timer = GLib.timeout_add_seconds(1, self._watchdog)
                    if self.pipeline.set_state(Gst.State.PLAYING) == Gst.StateChangeReturn.FAILURE:
                        raise RuntimeError("GStreamer не перешёл в PLAYING; проверьте плагины и доступ к NVDEC/NVENC.")
                    if not self.stopping:
                        self.loop.run()
                finally:
                    if self.pipeline is not None:
                        self.pipeline.set_state(Gst.State.NULL)
                        self.pipeline = None
                    if timer is not None:
                        GLib.source_remove(timer)
                    if bus is not None:
                        if handler is not None:
                            bus.disconnect(handler)
                        bus.remove_signal_watch()
                if self.error is not None:
                    raise self.error
                if self.stopping or not self.restart:
                    break
        finally:
            signal.signal(signal.SIGTERM, previous_sigterm)
