#include "config.h"
#include <stdlib.h>

/* Parse and validate a case without assembling a large model. */
int main(int argc, char **argv)
{
    PeecConfig config;
    config_set_defaults(&config);
    int result = config_parse_cli(&config, argc, argv);
    config_free(&config);
    return result < 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
