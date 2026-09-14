# Windows: только Fortran, Visual Studio и Intel ifx

## Установка

1. Установите Visual Studio 2022 с компонентом **Desktop development with C++**
   (Разработка классических приложений на C++): он предоставляет Windows SDK и
   компоновщик Microsoft. Исходники C этого проекта собирать не нужно.
2. Установите Intel Fortran Compiler **ifx** с интеграцией в установленную Visual Studio.
3. Установите Intel oneAPI Math Kernel Library (**oneMKL**): решатель вызывает LAPACK.
   При установке выбирайте совместимые версии компонентов. Если Visual Studio
   установлена позже компилятора, повторите установку/восстановление интеграции Intel.
4. Распакуйте ветку codex/fortran-port, например в C:\PEEC. В этой папке должны
   непосредственно находиться fortran, cases и meshes.

Python, GCC, WSL и исходный C-решатель для этой сборки и запуска не нужны.
oneMKL — внешняя численная библиотека; все исходники самого решателя имеют расширение .f90.

## Быстрая сборка

Откройте **Intel oneAPI command prompt for Visual Studio 2022** из меню Пуск:

```bat
cd /d C:\PEEC
fortran\windows\build_ifx.bat
fortran\windows\run.bat
```

Получится build\windows\PeecSolverFortran.exe. Скрипт собирает только модули
Fortran и main.f90; /Qopenmp включает OpenMP, /Qmkl:sequential подключает LAPACK
из последовательной oneMKL. Не включайте ILP64 или изменение стандартного INTEGER
на 8 байт: интерфейсы LAPACK рассчитаны на LP64.

Для прямого запуска из корня:

```bat
build\windows\PeecSolverFortran.exe
build\windows\PeecSolverFortran.exe --config cases/run_scattering.cfg
```

## Сборка в окне Visual Studio

1. File → New → Project → Intel Fortran → Console Application → Empty Project.
   Если шаблон отсутствует, проверьте установку интеграции Intel.
   Назовите проект PeecSolverFortran и создайте его, например, в C:\PEEC\vs.
2. Выберите платформу **x64** и компилятор Intel Fortran **ifx**.
3. Add → Existing Item: добавьте все десять файлов из fortran/src и
   fortran/app/main.f90. Не добавляйте fortran/test/test_core.f90:
   это отдельная тестовая программа со своим PROGRAM.
4. В свойствах проекта для нужной конфигурации (или All Configurations):
   - Fortran → Command Line → Additional Options: /Qopenmp /Qmkl:sequential
   - Debugging → Working Directory: C:\PEEC (ваш корень распакованного проекта)
   - Debugging → Command Arguments: оставить пустым для меню.
   Параметр oneMKL также можно выбрать через Fortran → Libraries →
   Intel oneAPI Math Kernel Library → Sequential.
5. Build → Build Solution. Запуск: Ctrl+F5; отладка: F5.
   Не меняйте размер стандартного INTEGER и правила именования внешних процедур.

Visual Studio определяет зависимости USE между модулями. Для ручной сборки порядок:
peec_base.f90 → peec_linalg.f90 → peec_quadrature.f90 → peec_mesh.f90 → peec_config.f90 → peec_model.f90 → peec_physics.f90 → peec_transient.f90 → peec_aperture_json.f90 → peec_io.f90 → main.f90.

## Выбор задачи в main.f90

В начале fortran/app/main.f90:

```fortran
integer, parameter :: selected_task=0
character(*), parameter :: custom_config='cases/run_scattering.cfg'
```

| selected_task | Задача / конфигурация |
|---|---|
| 0 | Меню в консоли; в меню 0 завершает программу |
| 1 | Информация о meshes/plate.msh |
| 2 | Рассеяние: cases/run_scattering.cfg |
| 3 | ЭПР: cases/run_rcs.cfg |
| 4 | Молния, замкнутое тело: fortran/cases/lightning_closed.cfg |
| 5 | Молния, два этапа: cases/run_lightning.cfg |
| 6 | Молния, slot-cells: cases/run_lightning_slot_cells.cfg |
| 7 | Файл custom_config; в меню путь вводится с клавиатуры |

Например, selected_task=3 запускает ЭПР после пересборки. Частоту, сетку,
шаг по времени и остальные параметры меняйте в соответствующем .cfg —
для изменения .cfg пересборка не нужна. Аргументы командной строки имеют
приоритет над selected_task. В меню за один запуск выполняется одна задача.

Начните с пункта 1. Пункты 2–6 используют более дорогие сетки body_*.
Для быстрого расчёта выберите пункт 7 и fortran/cases/plate_scattering.cfg.
Параметры щели и номера узлов в исходных конфигурациях требуют проверки
для конкретной физической модели.

## Пути и ошибки

- Запускайте из корня проекта; иначе cases/meshes не найдутся.
- В .cfg можно использовать /; пути с пробелами заключайте в кавычки.
  В пункте 7 меню вводите сам путь без окружающих кавычек.
- Выходные каталоги создаются автоматически на Windows и Linux.
  В Windows имена выходных путей не должны содержать символы "%!&|<>^.
- Ошибки LNK2019 с dgetrf/zgetrf и другими LAPACK-функциями:
  проверьте установку oneMKL и /Qmkl:sequential.
- Сообщение об отсутствующей Intel DLL: запускайте из oneAPI command prompt
  или установите соответствующие Intel runtime-библиотеки.
- Нет ifx в командной строке: используйте именно oneAPI command prompt.

## Проверка и документация

Численный перенос проверяется GNU Fortran в GitHub Actions; сборку в самой
Visual Studio/Intel ifx необходимо проверить на Windows. Не считайте инструкцию
свидетельством уже выполненной проверки Intel-компилятора.

Официальные инструкции Intel:
- [Выбор ifx в Visual Studio и командной строке](https://www.intel.com/content/www/us/en/developer/articles/training/intel-fortran-compiler-in-ms-visual-studio.html)
- [Параметр /Qmkl](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2025-0/qmkl-qmkl.html)
