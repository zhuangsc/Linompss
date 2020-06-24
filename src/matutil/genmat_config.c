#include "genmat_config.h"


#include "genmat.h"
#include "matfprint.h"
#include "strutil.h"


#if HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_LIBCONFIG
#include "libconfig.h"
#include <stdlib.h>
#include <string.h>
#endif


#ifdef SINGLE_PRECISION

const char *matmap[8] = {"spd", "covar", "rand", "eye", "diagdom", "ones", "zero", "costas"};

#define __genmat_config			sgenmat_config

#else 

extern const char *matmap[8];
#define __genmat_config			dgenmat_config

#endif


void* __genmat_config(char *app, char *id, int m, int n, int ldim, ompssmat_t deflt) 
{
	int dump = 0;
	double pars[4] = {2.0, 2.0, 2.0, 2.0};
	char tmpstr[128];

#if HAVE_LIBCONFIG
	const char *cfgfname = getenv("LINOMPSS_CFG");

	if ( cfgfname ) {
		config_t cfg;
		config_init(&cfg);

		snprintf(tmpstr, 128, "%spar", id);

		if ( config_read_file(&cfg, cfgfname)) {
			config_setting_t *setting = config_lookup(&cfg, app);
			if ( setting != NULL ) {
				char *matname = NULL;
				config_setting_lookup_string(setting, id, &matname);

				if ( matname != NULL ) {
					deflt = str2int(matname, matmap);
					printf("genmat: using %s %s (%i)\n", id, matname, deflt);
				}

				config_setting_t *paramsetting = config_setting_get_member(setting, tmpstr);
				if ( paramsetting ) {
			   		int count = config_setting_length(paramsetting);

					printf("genmat: parameters ");
					int i;
    				for( i = 0; i < count; ++i ) {
				  		double paramobj = config_setting_get_float_elem(paramsetting, i);
						pars[i] = paramobj;
						printf("%f ", paramobj);
					}
					printf("\n");
				}

				config_setting_t *dumpsetting = config_setting_get_member(setting, "dump");
				if ( dumpsetting ) {
			   		int count = config_setting_length(dumpsetting);

					int i;
    				for( i = 0; i < count; ++i ) {
				  		const char *dumpobj = config_setting_get_string_elem(dumpsetting, i);
						dump = ( strcmp(dumpobj, id) == 0 );  
						if ( dump ) break;
					}
				} 
			} else {
				printf("genmat: no settings for %s\n", app);
			}

			config_destroy(&cfg);
		} else {
			fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),config_error_line(&cfg), config_error_text(&cfg));
			config_destroy(&cfg);
//			fprintf(stderr, "genmat: cannot read config file %s\n", cfgfname); 
		}
	}
#endif

	void *res = NULL;

	switch(deflt) {
	case GENMATCONF_SPD:
		res = GENMAT_SPD(m, ldim);
		break;
	case GENMATCONF_COVAR:
		res = GENMAT_COVARIANT(m, pars[0], pars[1]); 
		break;
	case GENMATCONF_RAND:
		res = GENMAT(m, n, ldim); 
		break;
	case GENMATCONF_EYE:
		res = GENMAT_EYE(m, ldim); 
		break;
	case GENMATCONF_DIAGDOM:
		break;
	case GENMATCONF_ONES:
		break;
	case GENMATCONF_ZERO:
		res = GENMAT_ZERO(m, ldim); 
		break;
	case GENMATCONF_COSTAS:
		res = GENMAT_COSTAS(m, n); 
		break;
	default:
		break;
	}

	if ( dump ) {
		snprintf(tmpstr, 128, "%s.mm", id);
		FPRINT_DENSE2MM(tmpstr, id, m, n, res, ldim);
	}

	return res;
}
