
/**
 *  \file   ban-propag.c
 *  \brief  Deterministic BAN channel propagation layer
 *  \author Paul Ferrand, Javier Cuadrado
 *  \date   2010
 **/
#include <include/modelutils.h>
#include <math.h> // Use fmod(x,y)
#include "libcsv/csv.h"


/* ************************************************** */
/* ************************************************** */
model_t model =  {
    "BAN propagation model",
    "Paul Ferrand, Javier Cuadrado",
    "0.1",
    MODELTYPE_PROPAGATION, 
    {NULL, 0}
};

#define VERBOSE

#ifdef VERBOSE
#define VERB(X) X
#else
#define VERB(X)
#endif

#ifndef BAN_NODES_ACRONYMS
#define BAN_NODES_ACRONYMS
#define BN_POS_UNDEFINED -1
#define BN_POS_HIP 0
#define BN_POS_BACK 1
#define BN_POS_RTHIGH 2
#define BN_POS_RFOOT 3
#define BN_POS_LTHIGH 4
#define BN_POS_LFOOT 5
#define BN_POS_TORSO 6
#define BN_POS_RARM 7 
#define BN_POS_RHAND 8
#define BN_POS_LARM 9
#define BN_POS_LHAND 10
#define BN_POS_REAR 11
#define BN_POS_LEAR 12
#endif

#ifndef BAN_FADING_MODEL
#define BAN_FADING_MODEL
#define BF_NONE 	0
#define BF_RICE 	1
#define BF_NAKAGAMI	2
#define BF_RAYLEIGH	3
#endif

/* ************************************************** */
/* ************************************************** */
typedef struct samples 
{
	int number;
	int number_read;
	double time;
	double * values;
} t_samples;

struct entitydata {
	int fading_model;
	int node_count;
	t_samples * table;
	int table_init;
};

typedef struct csv_params
{
	int src;
	int dst;
	char * filename;
	int n_samples;
	double time;
	struct entitydata * entitydata; /* for the callback function, not very pretty... */
} t_csv_params;

/* ************************************************** */
/* *************** Normal distribution ************** */
/* ************************************************** */
double normal (double avg, double deviation) {
    return (avg + deviation * cos(2*M_PI*get_random_double()) * sqrt(-2.0 * log(get_random_double())));
}

/* ************************************************** */
/* ************ Manipulate sample tables ************ */
/* ************************************************** */

int init_tables(struct entitydata * entitydata)
{
	if (entitydata->table_init)
	{
		VERB(fprintf(stderr,"[ban-propag] Tables initialized. Use free_tables and this function afterwards to reset the values (init_tables).\n"));
		return -1;
	}
	
	VERB(fprintf(stderr,"[ban-propag] Initializing tables (init_tables).\n"));
	if (entitydata->node_count > 0)
	{
		entitydata->table = malloc( entitydata->node_count * entitydata->node_count * sizeof(t_samples) );
		entitydata->table_init = 1;
		return 0;
	}
	else
	{
		VERB(fprintf(stderr,"[ban-propag] node_count must be strictly positive (init_tables).\n"));
		return -1;
	}
}

__inline__
t_samples * get_samples_t( struct entitydata * entitydata, int src, int dest)
{
	return (entitydata->table + src * entitydata->node_count + dest);
}


int free_tables(struct entitydata * entitydata)
{
	int i;

	if (!entitydata->table_init)
	{
		VERB(fprintf(stderr,"[ban-propag] Tables not initialized (free_tables). Exiting ... \n"));
		return -1;
	}
	
	VERB(fprintf(stderr,"[ban-propag] Freeing tables (init_tables).\n"));
	if (entitydata->node_count > 0)
	{
		for (i = 0; i < entitydata->node_count * entitydata->node_count; i++)
		{
			free((entitydata->table + i)->values);			
		}
		free(entitydata->table);
		entitydata->table_init = 0;
		return 0;
	}
	else
	{
		VERB(fprintf(stderr,"[ban-propag] node_count must be strictly positive (free_tables).\n"));
		return -1;
	}
}

void reset_csv_params( t_csv_params * csv_params )
{
	csv_params->src = -1;
	csv_params->dst = -1;
	csv_params->filename = NULL;
	csv_params->n_samples = -1;
	csv_params->time = -1.0;	
}

/* ************************************************** */
/* ************** Callbacks for libcsv ************** */
/* ************************************************** */

/* 
	"cb1" for sample data
	*p should contain a pointer to the appropriate samples structure (t_samples *)
	p->number_read is incremented upon successful read and convert operation.
*/
void _links_field_read(void *s, size_t i, void *p) 
{
	char buf[64];
	strncpy( buf, (char *)s, i );
	buf[i] = '\0';
	
	t_samples * samples = (t_samples *)p;
	if (samples->number_read == samples->number)
	{
		VERB(fprintf(stderr,"[ban-propag] Too many samples (%d announced in the configuration file) (_links_field_read). \n", samples->number));
		VERB(fprintf(stderr,"[ban-propag] To prevent an overflow the sample is discarded (_links_field_read). \n"));
		return;
	}
	samples->values[samples->number_read++] = atof( buf );
}	

/* 
	"cb2" for sample data
	*p should contain a pointer to the appropriate samples structure (t_samples *)
*/
void _links_entry_read(int i, void *p)
{
	VERB(t_samples * samples = (t_samples *)p);
	VERB(fprintf(stdout,"[ban-propag] Link data read (_links_entry_read).  \n"));
	VERB(fprintf(stdout,"[ban-propag] Entries : %d.  \n", samples->number_read ));
}

/* 
	"cb1" for the configuration file
	*p should contain a pointer to the appropriate param structure (t_csv_params *)
*/
void _config_field_read(void *s, size_t i, void *p) 
{
	/* Temporary char buffer, copying the seemingly not null-terminated string
		from the library before conversion */
	char buf[256];
	strncpy( buf, (char *)s, i );
	buf[i] = '\0';
	
	/*VERB(fprintf(stdout,"Param read : %s\n",buf));*/
	
	t_csv_params * csv_params = (t_csv_params *)p;
	if (csv_params->src == -1) /* First param (source) not set, we consider we read the first parameter */
	{
		csv_params->src = atoi( buf );
		return;
	}
	
	if (csv_params->dst == -1) /* Second param (destination)*/
	{
		csv_params->dst = atoi( buf );
		return;
	}
	if (csv_params->filename == NULL) /* Third param : filename */
	{
		if (csv_params->filename)
			free(csv_params->filename);
		csv_params->filename = malloc( i * sizeof(char));
		strncpy( csv_params->filename, buf, i );
		csv_params->filename[i] = '\0';
		return;
	}
	if (csv_params->n_samples == -1) /* Fourth param : number of samples */
	{
		csv_params->n_samples = atoi( buf );
		return;
	}
	if (csv_params->time == -1.0) /* Fifth param : time sampled */
	{
		csv_params->time = atof( buf );
		return;
	}
}

/* 
	"cb2" for the configuration file
	*p should contain a pointer to the appropriate param structure (t_csv_params *)
*/
void _config_entry_read(int i, void *p)
{
	t_csv_params * csv_params = (t_csv_params *)p;
	t_samples * samples;
	FILE * fp;
	/* char buf; */
	char buf[1024]; 
	size_t bytes_read;
	struct csv_parser parser;
	
	VERB(fprintf(stdout,"[ban-propag] Configuration data read (_config_entry_read).  \n"));
	VERB(fprintf(stdout,"[ban-propag] Summary for the link : \n"));
	
	if (csv_params->src == -1) 
	{
		VERB(fprintf(stderr,"[ban-propag] Source not set! \n"));
		goto error;
	}
	else
	{
		VERB(fprintf(stdout,"[ban-propag] Source : %d \n", csv_params->src));
	}	
	if (csv_params->dst == -1) 
	{
		VERB(fprintf(stderr,"[ban-propag] Destination not set! \n"));
		goto error;
	}
	else
	{
		VERB(fprintf(stdout,"[ban-propag] Destination : %d \n", csv_params->dst));
	}	
	if (!csv_params->filename) 
	{
		VERB(fprintf(stderr,"[ban-propag] File name not set! \n"));
		goto error;
	}
	else
	{
		VERB(fprintf(stdout,"[ban-propag] File name : %s \n", csv_params->filename));
	}
	if (csv_params->n_samples == -1)
	{
		VERB(fprintf(stderr,"[ban-propag] Sample count not set! \n"));
		goto error;
	}
	else
	{	
		VERB(fprintf(stdout,"[ban-propag] Sample count : %d \n", csv_params->n_samples));
	}
	if (csv_params->time == -1.0f) 
	{
		VERB(fprintf(stderr,"[ban-propag] Sampled time not set! \n"));
		goto error;
	}
	else
	{
		VERB(fprintf(stdout,"[ban-propag] Sampled time : %f \n", csv_params->time));
	}
	
	samples = get_samples_t(csv_params->entitydata, csv_params->src, csv_params->dst);
	if (!samples->values)
		samples->values = malloc( csv_params->n_samples * sizeof(double) );
	samples->time = csv_params->time;
	samples->number = csv_params->n_samples;
	
	if( (fp = fopen(csv_params->filename, "rb")) )
	{
		csv_init(&parser, 0);
		while( ( bytes_read = fread(buf, 1, 1024, fp) ) > 0 ) 
		/* while( ( buf=fgetc(fp) ) != EOF) */
		{
			if (csv_parse(&parser, buf, bytes_read, _links_field_read, _links_entry_read, samples) != bytes_read) 
			{
				VERB(fprintf(stderr, " [ban-propag] Error: %s (_config_entry_read)\n", csv_strerror(csv_error(&parser))));
				goto error;
			}
		}
		csv_fini(&parser, _links_field_read, _links_entry_read, samples);
		csv_free(&parser);
		fclose(fp);
	}
	else
	{
		VERB(fprintf(stderr,"[ban-propag] Unable to open the sample file (%s).\n", csv_params->filename));
		VERB(fprintf(stderr,"[ban-propag] The file may be inexistent, or has too many read handles already opened.\n"));
		goto error;
	}	
	reset_csv_params(csv_params);
	return;
	
	error:
	if (csv_params->filename)
			free(csv_params->filename);
	csv_free(&parser);
	reset_csv_params(csv_params);
	VERB(fprintf(stderr,"[ban-propag] An error occured, ignoring this record (_config_entry_read). \n"));
}

/* ************************************************** */
/* ************************************************** */
int init(call_t *c, void *params) 
{
	struct entitydata * entitydata = malloc(sizeof(struct entitydata));
	param_t * param;
	FILE * fp;
	struct csv_parser parser;
	t_csv_params csv_params;
	/* char buf; */
	char buf[1024]; 
	size_t bytes_read;
/* Not using a 1024 char buffer because he seems to be filled with junk from time to time o_o */
	
	
	/* default values */
	entitydata->table_init = 0;
	entitydata->fading_model = BF_NONE;
	entitydata->node_count = get_node_count();
	reset_csv_params(&csv_params);
	/* get parameters */
	das_init_traverse(params);
	while ((param = (param_t *) das_traverse(params)) != NULL) {
        if (!strcmp(param->key, "data_description_file")) {
			fp = fopen( param->value, "rb" );
			if( fp ) 
			{
				VERB(fprintf(stdout,"[ban-propag] Reading the data description file (%s).\n", param->value));
				init_tables(entitydata);
				csv_params.entitydata = entitydata; /* pointer to the table for the callback functions */
				
				/* Read the csv data description file */
				/* for each line fill the tables */
				csv_init(&parser, 0);
				while( ( bytes_read = fread(buf, 1, 1024, fp) ) > 0 ) 
				/*while( ( buf=fgetc(fp) ) != EOF) */
				{
	 				if (csv_parse(&parser, buf, bytes_read, _config_field_read, _config_entry_read, &csv_params) != bytes_read) 
					{
						VERB( fprintf(stderr, "Error: %s\n", csv_strerror(csv_error(&parser))) );
      					goto error;
    				}
  				}
  				csv_fini(&parser, _config_field_read, _config_entry_read, &csv_params);
				csv_free(&parser);
				fclose(fp);
 			} 
			else 
			{
				VERB(fprintf(stderr,"[ban-propag] Unable to open the data description file (%s).\n", param->value));
				VERB(fprintf(stderr,"[ban-propag] The file may be inexistent, or has too many read handles already opened.\n"));
			}
        }
        if (!strcmp(param->key, "fading_model")) 
		{
			if (!strcmp(param->value, "rayleigh"))
			{
				entitydata->fading_model = BF_RAYLEIGH;
				VERB(fprintf(stderr,"[ban-propag] Fading set to Rayleigh.\n"));
			}
			else if (!strcmp(param->value, "rician"))
			{
				entitydata->fading_model = BF_RICE;
				VERB(fprintf(stderr,"[ban-propag] Fading set to Rician K.\n"));
			}
			else if (!strcmp(param->value, "nakagami"))
			{
				entitydata->fading_model = BF_NAKAGAMI;
				VERB(fprintf(stderr,"[ban-propag] Fading set to Nakagami-m.\n"));
			}
			else if (!strcmp(param->value, "none"))
			{
				entitydata->fading_model = BF_NONE;
				VERB(fprintf(stderr,"[ban-propag] Fading set to None.\n"));
			}
            else 
			{
				VERB(fprintf(stderr,"[ban-propag] Unknown fading type : (%s) !\n", param->value));
				goto error;
			}
        }
    }
    
    set_entity_private_data(c, entitydata);
    VERB(fprintf(stderr,"[ban-propag] Entity data set!\n"));
	return 0;

 error:
	 VERB(fprintf(stderr,"[ban-propag] Something went wrong! Aborting the initialization...\n"));
	if (entitydata->table_init)
		free_tables(entitydata);
	free(entitydata);
    return -1;
}

int destroy(call_t *c) {
	VERB(fprintf(stderr,"[ban-propag] Destroying the propagation framework!\n"));
	struct entitydata * edata = get_entity_private_data(c);
	if (edata->table_init)
		free_tables(edata);
    free(edata);
    return 0;
}



/* ************************************************** */
/* ************************************************** */
double compute_ban_fading(double pathloss, int fading_model)
{
    double fading_env = 0.0;
    double K_db = 0.0;
    double K = 0.0;
    double m = 1.0;
    double sigma = 1/sqrt(2);
    switch(fading_model)
    {
		case BF_NONE: 
			fading_env = 0.0; 
			break;
		case BF_RAYLEIGH:
			fading_env = hypot(normal(0,sigma),normal(0,sigma)); 
		case BF_RICE: 
			K_db = 0.43*pathloss + 6*(get_random_double()-0.5);
			K = pow(10.0,(K_db/10.0));
			fading_env =  hypot(sqrt(K/(K+1)) * sigma + sqrt(1/(K+1)) * 
								normal(0,sigma),
						sqrt(1/(K+1))*normal(0,sigma)); 
			break;
		case BF_NAKAGAMI: 
			fading_env = 0.0; 
			break;
		default : fading_env = 0.0; 
    }
    VERB(fprintf(stdout,"[ban-propag] Fading enveloppe r.v. : %f.\n", fading_env));
    if (fading_env > 0.0)
	    return 20*log10(fading_env);
	else
		return 0.0;
}

/* ************************************************** */
/* ************************************************** */
double propagation(call_t *c, packet_t *packet, nodeid_t src, nodeid_t dst, double rxdBm) 
{
    struct entitydata * entitydata = get_entity_private_data(c);
    t_samples * samples;
    double sim_time;
    double pos;
    int index;
    double pathloss, fading;
    double rx_dBm = rxdBm;
    
    if ( !(samples = get_samples_t( entitydata, src, dst)) )
    {
    	VERB(fprintf(stdout,"[ban-propag] Unallocated samples structure. Aborting (propagation).\n"));	
    	return rx_dBm;
    }
    
#define NANO (1000000000);
	sim_time = get_time()/NANO;
	fprintf(stdout,"[ban-propag] Time called : %f s (propagation).\n", sim_time);
    
    pos = fmod(sim_time, samples->time);
    index = (pos/samples->time) * samples->number;
    
    if ( (index < 0) || (index >= samples->number) )
    {
    	VERB(fprintf(stdout,"[ban-propag] Sample index out of bounds!.\n"));
    	VERB(fprintf(stdout,"[ban-propag] Simulation time : %f.\n", sim_time));
    	VERB(fprintf(stdout,"[ban-propag] Modulo sample time : %f.\n", pos));
    	VERB(fprintf(stdout,"[ban-propag] Index/samples : %d/%d.\n", index, samples->number));
    	return rx_dBm;    	    	
    }
    
    pathloss = samples->values[index];
    rx_dBm -= pathloss;
    fading = compute_ban_fading(pathloss, entitydata->fading_model);
    rx_dBm += fading;    
    
    VERB(fprintf(stdout,"[ban-propag] P_tx (dBm) : %f.\n", rxdBm));
    VERB(fprintf(stdout,"[ban-propag] Pathloss (shadowing) (dB) : %f.\n", pathloss));
    VERB(fprintf(stdout,"[ban-propag] Fading (dB) : %f.\n", fading));
    VERB(fprintf(stdout,"[ban-propag] P_rx (dBm) : %f.\n", rx_dBm));
    return rx_dBm;
}


/* ************************************************** */
/* ************************************************** */
propagation_methods_t methods = {propagation};
