#include "ompss_apps.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_LIBCONFIG
#include "libconfig.h"
#include <stdlib.h>
#include <string.h>
#include "strutil.h"
#endif


const char *ompssapp_map[6] = {"CG", "CGMOD1", "CGPROF", "POTRF", "LU", "QR"};


const char *ompssapp_CG = "CG";
const char *ompssapp_CGPROF = "CGPROF";
const char *ompssapp_CGMOD1 = "CGMOD1";
const char *ompssapp_POTRF = "POTRF";
const char *ompssapp_LU = "LU";
const char *ompssapp_QR = "QR";
const char *ompssapp_ITREF = "ITREF";


int ompssapp_readid(char *app, char *id, ompssapp_t deflt) 
{
#if HAVE_LIBCONFIG
	const char *cfgfname = getenv("LINOMPSS_CFG");

	if ( cfgfname ) {
		config_t cfg;
		config_init(&cfg);

		if ( config_read_file(&cfg, cfgfname)) {
			config_setting_t *setting = config_lookup(&cfg, app);
			if ( setting != NULL ) {
				char *idname = NULL;
				config_setting_lookup_string(setting, id, &idname);

				if ( idname != NULL ) {
					deflt = str2int(idname, ompssapp_map);
					printf("ompssapp: using %s %s (%i)\n", id, idname, deflt);
				}
			}

			config_destroy(&cfg);
		} else {
			fprintf(stderr, "ompssapp: cannot read config file %s\n", cfgfname); 
		}
	}
#endif

	return deflt;
}
