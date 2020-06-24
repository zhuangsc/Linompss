#include "ompss_options.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_LIBCONFIG
#include "libconfig.h"
#include <stdlib.h>
#include <string.h>
#include "strutil.h"
#endif


const char *ompssopt_map[6] = {"no", "yes"};


const char *ompssopt_NO = "no";
const char *ompssopt_YES = "yes";


int ompssopt_read(char *app, char *opt, ompssopt_t deflt) 
{
#if HAVE_LIBCONFIG
	const char *cfgfname = getenv("LINOMPSS_CFG");

	if ( cfgfname ) {
		config_t cfg;
		config_init(&cfg);

		if ( config_read_file(&cfg, cfgfname)) {
			config_setting_t *setting = config_lookup(&cfg, app);
			if ( setting != NULL ) {
				char *optname = NULL;
				config_setting_lookup_string(setting, opt, &optname);

				if ( optname != NULL ) {
					deflt = str2int(optname, ompssopt_map);
					printf("ompssopt: using %s %s (%i)\n", opt, optname, deflt);
				}
			}

			config_destroy(&cfg);
		} else {
			fprintf(stderr, "ompssopt: cannot read config file %s\n", cfgfname); 
		}
	}
#endif

	return deflt;
}

