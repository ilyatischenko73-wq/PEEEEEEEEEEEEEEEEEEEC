/* This narrow probe calls the actual static file parser from the driver.
 * When that parser becomes public, replace this include with its header.
 */
#define main peec_driver_main
#include "../src/main.c"
#undef main

/* Test address lifetime and undefined behavior here, not ownership cleanup.
 * LeakSanitizer is disabled because this probe does not run driver teardown;
 * it also makes the probe usable in containers without /proc access.
 */
const char *__asan_default_options(void) { return "detect_leaks=0"; }
int __lsan_is_turned_off(void) { return 1; }
const char *__ubsan_default_options(void) { return "halt_on_error=1"; }

static int matches_live_string(const char *value, const char *expected)
{
    if (value == NULL) return 0;
    /* Force actual, instrumented reads, including the terminating byte. */
    const volatile char *live = value;
    for (size_t i = 0; ; ++i) {
        if (live[i] != expected[i]) return 0;
        if (expected[i] == '\0') return 1;
    }
}

int main(void)
{
    PeecConfig config;
    config_set_defaults(&config);
    if (config_parse_file_as_cli(
            &config, "test_config_lifetime", "tests/fixtures/file_mode.cfg") != 0) {
        fprintf(stderr, "File configuration must parse successfully.\n");
        return EXIT_FAILURE;
    }
    int ok = matches_live_string(config.mesh_file, "meshes/plate.msh")
        && matches_live_string(config.matrix_directory, "tests-out/matrices")
        && matches_live_string(config.vtk_directory, "tests-out/vtk")
        && matches_live_string(config.output_file, "tests-out/rcs.csv");
    puts(ok ? "PASS config_file_lifetime" : "FAIL config_file_lifetime");
    config_free(&config);
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
