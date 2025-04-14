#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<glib.h>
#include<regex.h>
#include<limits.h>
#include<stddef.h>
#include<unistd.h>
#include<ctype.h>
#include<stdbool.h>
#include<assert.h>

#define MAX_CHROMOSOMES 23
#define TOTAL_MARKERS 500
#define BUFFER_SIZE 60000
#define MAX_HAPSEQS 1024
#define MAX_INDIVIDUALS 11000
#define MAX_FILENAME 100
#define MAX_INDIVIDUAL_NAME 10
#define NO_MODELS 6
#define add_hap 1
#define remove_hap 0
#define model_a 1
#define model_b 2
#define model_c 3
#define model_d 4
#define model_e 5
#define model_f 6

#undef DEBUG

#define PLL_STRING(x) #x
#define PLL_C2S(x) PLL_STRING(x)

#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define VERSION_PATCH 0

#define PROG_VERSION "v" PLL_C2S(VERSION_MAJOR) "." PLL_C2S(VERSION_MINOR) "." \
        PLL_C2S(VERSION_PATCH)

#define RESET "\x1B[0m"
#define BOLD "\x1B[1m"      //Bold Text Formula...


/* const int PROG_BAR_LENGTH = 100; */
/* void printProgress(double percentage) { */
/*     int val = (int) (percentage * 100); */
/*     int lpad = (int) (percentage * PBWIDTH); */
/*     int rpad = PBWIDTH - lpad; */
/*     printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, ""); */
/*     fflush(stdout); */
/* }    */


GError *error = NULL;
GIOChannel *infile_chrom = NULL;
GIOChannel *infile_A = NULL;
GIOChannel *infile_B = NULL;
GIOChannel *infile_indv = NULL;
GIOChannel *outfile = NULL;

char prog_name[] = "mongrail2";
char version[] = "mongrail2";
char progheader[100];

void fillheader(const char version[])
{
  snprintf(progheader, 80,
           "%s %s %s",
           prog_name, PROG_VERSION, version);
}

void show_header()
{
  fprintf(stdout, "\n%s\n", progheader);
  fprintf(stdout, "https://github.com/Mongrail-2-0/\n");
  fprintf(stdout,"\n");
}

static void program_help()
{
  printf("Usage: ./mongrail2 -c FILENAME -A FILENAME -B FILENAME -i FILENAME [-o FILENAME]\n\n");
  printf("%s   Required Arguments: %s\n"
	 "\t-c, --chromosome=FILE\n\t\tTake chromosome information from FILE.\n\n"
	 "\t-A, --Apop=FILE\n\t\tTake population A haplotype information from FILE.\n\n"
	 "\t-B, --Bpop=FILE\n\t\tTake population B haplotype information from FILE.\n\n"
	 "\t-i, --input=FILE\n\t\tTake sampled haplotypes from FILE.\n\n"
	 "%s   Optional arguments: %s\n"
	 "\t-o, --output=FILE\n\t\tPrint output to FILE. If this option is not used the default output filename is 'out.txt'.\n\n"
	 "\t-V, --verbose\n\t\tPrint verbose output to screen.\n\n", BOLD, RESET, BOLD, RESET
	 );
}

static void print_msg()
{
    printf("Usage: ./mongrail2 -c FILENAME -A FILENAME -B FILENAME -i FILENAME [-o FILENAME]\n Try './mongrail2 -h' for more information.\n\n");
  /* printf("Usage: mongrail2 [OPTION]... FILE\n Try 'mongrail2 -h' for more information.\n\n"); */

}


void update_bar(int percent_done, unsigned int chrom_index)
{
  /* int num_chars = percent_done * PROG_BAR_LENGTH / 100; */
  /* printf("["); */
  /* for (int i=0; i < num_chars; i++) */
  /*   { */
  /*     printf("\r#"); */
  /*   } */
  /* for (int i=0; i < PROG_BAR_LENGTH-num_chars; i++) */
  /*   { */
  /*     printf(" "); */
  /*   } */
  /* printf("] %d%% Done", percent_done); */
  printf("\rCalculating the log-likelihood for chromosome %u: %d%% Done", chrom_index, percent_done);
  fflush(stdout);
}


typedef struct chromosome_all_info
{
  unsigned int chrom_id;
  double chrom_length;
  double chrom_recom_rate;
  unsigned int n_loci;
  double *markers;
  GHashTable* hashA;
  GHashTable* hashB;
  unsigned int *haplotype_1;
  unsigned int *haplotype_2;
  struct chromosome_all_info *next;
} chrom_data;


chrom_data* make_data_node(unsigned int index, double length, double recom_rate, unsigned int no_loci, double *loci)
{
  chrom_data *new_data_node;
  if((new_data_node = malloc(sizeof(chrom_data))) == NULL)
    { fprintf(stderr, "Oops, out of memory!"); exit(1);}

  new_data_node->chrom_id = index;
  new_data_node->chrom_length = length;
  new_data_node->chrom_recom_rate = recom_rate;
  new_data_node->n_loci = no_loci;
  new_data_node->markers = malloc(no_loci * sizeof(double));
  assert(new_data_node->markers != NULL);
  for(int i = 0; i < no_loci; i++)
    {
      new_data_node->markers[i] = loci[i];
    }
  new_data_node->hashA = NULL;
  new_data_node->hashB = NULL;
  new_data_node->haplotype_1 = NULL;
  new_data_node->haplotype_2 = NULL;
  new_data_node->next = NULL;
  return(new_data_node);
}

void string_to_markers(char *string, unsigned int no_loci, double **markers)
{
  char *second_regexString = "([0-9]+)";
  char *second_string;                                                
  second_string = string; 
  regex_t *second_regexCompiled = (regex_t*)malloc(sizeof(regex_t));
  assert(second_regexCompiled != NULL);
  regmatch_t *pmatch;                                                 
  if(regcomp(second_regexCompiled, second_regexString, REG_EXTENDED|REG_NEWLINE))             
    {                                                                 
      printf("Could not compile regular expression.\n");              
      exit(1);                                                        
    }                                                                 
  pmatch = (regmatch_t*)malloc(sizeof(regmatch_t)*second_regexCompiled->re_nsub);
  assert(pmatch != NULL);
  if(regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0))  
    {                                                                 
      printf("Nothing matched with ""%s""\n", second_string);         
      exit(1);                                                        
    }   
  unsigned int n_loci = 0;                                            
  double *loci_position;                                              
  loci_position = malloc(TOTAL_MARKERS * sizeof(double));
  assert(loci_position != NULL);
  do {                                                                
    if (pmatch[0].rm_so != -1) {        /* The regex is matching part\
					   of a string */                                                       
      char *submatch;                                                 
      double val;                                                     
      size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;            
      submatch = (char*)malloc(matchlen+1);
      assert(submatch != NULL);
      strncpy(submatch, second_string + pmatch[0].rm_so, matchlen+1); 
      submatch[matchlen]='\0';                                        
      val = atof(submatch);                                           
      loci_position[n_loci] = val;                                    
      free(submatch);                                                 
    };                                                                
    second_string += pmatch[0].rm_eo;   /* Restart from last match */ 
    n_loci++;                                                         
  } while(!regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch,0));
  *markers = malloc(n_loci * sizeof(double));
  assert(markers != NULL);
  for(int i = 0; i < n_loci; i++ )                                    
    {
      *(*markers + i) = loci_position[i];
#ifdef DEBUG
      printf("Locus[%d]: %lf\t %lf\n", i+1, loci_position[i],*(*markers + i));
#endif
    }  
  regfree(second_regexCompiled);                                      
  free(second_regexCompiled);                                         
  free(pmatch);
  free(loci_position);

}

void read_chrom(char *filename, GIOChannel *infile, unsigned int **chromosome_index, unsigned int **chromosome_no_loci, double **chromosome_length, double **chromosome_recom_rate, unsigned int *no_chromosomes, double **all_markers, unsigned int *no_total_markers)
{
  char *complete_file;
  if(g_io_channel_read_to_end(infile,&complete_file,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *source;
  source = complete_file;
  char * regexString = "([0-9]+)[[:blank:]]+([0-9]+\\.?[0-9]*)[[:blank:]]+([0-9\
]+\\.?[0-9]*)[[:blank:]]+([0-9]+)[[:blank:]]+(.+)";
  size_t maxGroups = 6;
                                                                                
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];

  unsigned int *chrom_index, *chrom_no_loci;
  double *chrom_length, *chrom_recom_rate;
  double *chrom_markers;
  chrom_index = calloc(MAX_CHROMOSOMES, sizeof(unsigned int));
  chrom_no_loci = calloc(MAX_CHROMOSOMES, sizeof(unsigned int));
  chrom_length = calloc(MAX_CHROMOSOMES, sizeof(double));
  chrom_recom_rate = calloc(MAX_CHROMOSOMES, sizeof(double));
  chrom_markers = calloc(MAX_CHROMOSOMES * TOTAL_MARKERS, sizeof(double));
  assert(chrom_index != NULL);
  assert(chrom_no_loci != NULL);
  assert(chrom_length != NULL);
  assert(chrom_recom_rate != NULL);
  assert(chrom_markers != NULL);


  if (regcomp(&regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    };

  /* printf("Reading chromosome information from file: %s ...\n", filename); /\* reading file status printout *\/ */
  
  unsigned int c = 0;
  unsigned int chrom_no = 0;
  unsigned int val_1, val_4;
  double val_2, val_3;
  int marker_index = 0;
  for(c = 0; c < MAX_CHROMOSOMES; c++)
    {
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
	{
	  unsigned int g = 0;
	  unsigned int offset = 0;
	  char sourceCopy[strlen(source) + 1];
	  strcpy(sourceCopy, source);
	  sourceCopy[groupArray[g].rm_eo] = 0;
	  for (g = 0; g < maxGroups; g++)
	    {
	      if (groupArray[g].rm_so == (size_t)-1)
		break;  // No more groups

	      if(g == 0)
		{
		  offset = groupArray[g].rm_eo;
		}
	      if(g == 1)
		{
		  char *submatch_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
		  val_1 = atoi(submatch_1);
		  chrom_index[chrom_no] = val_1;
#ifdef DEBUG
		  printf("Chromosome ID: %u\n\n", val_1);
#endif
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  char *submatch_2;
		  size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_2 = (char*)malloc(matchlen_2+1);
		  assert(submatch_2 != NULL);
		  strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		  submatch_2[matchlen_2]='\0';
		  val_2 = atof(submatch_2);
		  chrom_length[chrom_no] = val_2;
#ifdef DEBUG
		  printf("Chromosome length (in Mb): %lf\n\n", val_2);
#endif
		  free(submatch_2);
		}
	      if(g == 3)
		{
		  char *submatch_3;
		  size_t matchlen_3 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_3 = (char*)malloc(matchlen_3+1);
		  assert(submatch_3 != NULL);
		  strncpy(submatch_3, sourceCopy + groupArray[g].rm_so, matchlen_3+1);
		  submatch_3[matchlen_3]='\0';
		  val_3 = atof(submatch_3);
		  chrom_recom_rate[chrom_no] = val_3;
#ifdef DEBUG
		  printf("Recombination rate (in cM/Mb): %lf\n\n", val_3);
#endif
		  free(submatch_3);
		}
	      if(g == 4)
		{
		  char *submatch_4;
		  size_t matchlen_4 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_4 = (char*)malloc(matchlen_4+1);
		  assert(submatch_4 != NULL);
		  strncpy(submatch_4, sourceCopy + groupArray[g].rm_so, matchlen_4+1);
		  submatch_4[matchlen_4]='\0';
		  val_4 = atoi(submatch_4);
		  chrom_no_loci[chrom_no] = val_4;
#ifdef DEBUG
		  printf("Number of markers: %u\n\n", val_4);
#endif
		  free(submatch_4);
		}
	      if(g == 5)
	      	{
		  char *submatch_5;
		  size_t matchlen_5 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_5 = (char*)malloc(matchlen_5+1);
		  assert(submatch_5 != NULL);
		  strncpy(submatch_5, sourceCopy + groupArray[g].rm_so, matchlen_5+1);
		  submatch_5[matchlen_5]='\0';
		  /* printf("%s,\n\n",submatch_5); */
		  double *loci;
		  string_to_markers(submatch_5, chrom_no_loci[chrom_no], &loci);
		  for(int i = 0; i < val_4; i++)
		    {
		      chrom_markers[marker_index + i] = loci[i];
		    }
		  free(loci);
		  free(submatch_5);
	      	}
	    }
	  marker_index = marker_index + val_4;	  
	  source = source + offset;
	  chrom_no++;
	}
    }
  /* printf("# chrom:%u\n\n\n", chrom_no); */
  
  *chromosome_index = malloc(chrom_no * sizeof(unsigned int));
  *chromosome_no_loci = malloc(chrom_no * sizeof(unsigned int));
  *chromosome_length = malloc(chrom_no * sizeof(double));
  *chromosome_recom_rate = malloc(chrom_no * sizeof(double));
  *all_markers = malloc(marker_index * sizeof(double));
  assert(chromosome_index != NULL);
  assert(chromosome_no_loci != NULL);
  assert(chromosome_length != NULL);
  assert(chromosome_recom_rate != NULL);
  assert(all_markers != NULL);

  for(int i = 0; i < chrom_no; i++)
    {
      *(*chromosome_index + i) = chrom_index[i];
      *(*chromosome_no_loci + i) = chrom_no_loci[i];
      *(*chromosome_length + i) = chrom_length[i];
      *(*chromosome_recom_rate + i) = chrom_recom_rate[i];
    }
  for(int i = 0; i < marker_index; i++)
    {
      *(*all_markers + i) = chrom_markers[i];
    }

  *no_chromosomes = chrom_no;
  /* printf("Number of chromosomes read from file: %u\n", chrom_no); */
  *no_total_markers = marker_index;
  regfree(&regexCompiled);
  g_free(complete_file);
  free(chrom_index);
  free(chrom_no_loci);
  free(chrom_length);
  free(chrom_recom_rate);
  free(chrom_markers);
}

chrom_data* data_linked_list_creation(unsigned int *chromosome_index, unsigned int *chromosome_no_loci, double *chromosome_length, double *chromosome_recom_rate, unsigned int no_chromosomes, double *chromosome_markers)
{
  chrom_data* head, *current;
  int first = 0;
  int marker_index = 0;
  /* if(verbose) */
  /*   { */
  /*     printf("chrom index | #markers\n"); */
  /*     printf("%.*s\n",22,"---------------------------------------------------------------------"); */
  /*   } */

  for(int i = 0; i < no_chromosomes; i++)
    {
      double *loci;
      loci = chromosome_markers + marker_index;
      chrom_data *next_data_node = make_data_node(chromosome_index[i], chromosome_length[i], chromosome_recom_rate[i], chromosome_no_loci[i], loci);
      marker_index = marker_index + chromosome_no_loci[i];
      if(first == 0)
	{
	  head = next_data_node;
	  current = head;
	  /* if(verbose) */
	  /*   { */
	  /*     printf("%11u | %8u\n",current->chrom_id,current->n_loci); */
	  /*   } */
	  first = 1;
	}
      else
	{
	  current->next = next_data_node;
	  current = current->next;
	  /* if(verbose) */
	  /*   { */
	  /*     printf("%11u | %8u\n",current->chrom_id,current->n_loci); */
	  /*   } */
	}
    }
  /* printf("Finished reading the chromosome information file successfully!\n\n"); */
  return(head);
}


void hap_info_extract(char* string, unsigned int **sequence, unsigned int **count, unsigned int *no_sequence)
{
  char *second_string;
  second_string = string;
  char *second_regexString = "([0-1]+:[0-9]+)";
  regex_t *second_regexCompiled = (regex_t*)malloc(sizeof(regex_t));
  assert(second_regexCompiled != NULL);
  regmatch_t *pmatch;
  if(regcomp(second_regexCompiled, second_regexString, REG_EXTENDED|REG_NEWLINE))
    {                                                                            
      printf("Could not compile regular expression.\n");
      exit(1);
    }                                                                            
  pmatch = (regmatch_t*)malloc(sizeof(regmatch_t)*second_regexCompiled->re_nsub);
  assert(pmatch != NULL);
  if(regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0))
    {                                                                                       
      printf("Nothing matched with ""%s""\n", second_string);
      exit(1);
    }
  unsigned int n_hapseq = 0;
  unsigned int *hap_seq;
  unsigned int *hap_count;
  hap_seq = malloc(MAX_HAPSEQS * sizeof(unsigned int));
  hap_count = malloc(MAX_HAPSEQS * sizeof(unsigned int));
  assert(hap_seq != NULL);
  assert(hap_count != NULL);

  do {                                                                          
    if (pmatch[0].rm_so != -1) {        /* The regex is matching part of a string */
      char *submatch;
      size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;
      submatch = (char*)malloc(matchlen+1);
      assert(submatch != NULL);
      strncpy(submatch, second_string + pmatch[0].rm_so, matchlen+1);
      submatch[matchlen]='\0';
      /* printf("%s\n",submatch); */
      char *source;
      source = submatch;
      char * regexString = "([0-1]+):([0-9]+)";
      size_t maxGroups = 3;
      regex_t regexCompiled;
      regmatch_t groupArray[maxGroups];
      if (regcomp(&regexCompiled, regexString, REG_EXTENDED))                                      
        {                                                                                          
          printf("Could not compile regular expression.\n");
          exit(1);
        };
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
        {
          unsigned int g = 0;
          for (g = 0; g < maxGroups; g++)                                                          
            {                                                                                      
              if (groupArray[g].rm_so == (size_t)-1)                                               
                break;  // No more groups
              char sourceCopy[strlen(source) + 1];
              strcpy(sourceCopy, source);
              sourceCopy[groupArray[g].rm_eo] = 0;
	      if(g == 1)
		{
		  char *submatch_1;
		  long val_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
		  val_1 = strtol(submatch_1, NULL, 2);
		  hap_seq[n_hapseq] = (unsigned int) val_1;
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  char *submatch_2;
		  unsigned int val_2;
		  size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_2 = (char*)malloc(matchlen_2+1);
		  assert(submatch_2 != NULL);
		  strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		  submatch_2[matchlen_2]='\0';
		  val_2 = atoi(submatch_2);
		  hap_count[n_hapseq] = (unsigned int) val_2;
		  free(submatch_2);
		}
            }
        }
      regfree(&regexCompiled);
      free(submatch);
    };
    second_string += pmatch[0].rm_eo;   /* Restart from last match */
    n_hapseq++;
  } while(!regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0)); 
  *sequence = malloc(n_hapseq * sizeof(unsigned int));
  *count = malloc(n_hapseq * sizeof(unsigned int));
  assert(sequence != NULL);
  assert(count != NULL);

  for(int i = 0; i < n_hapseq; i++)
    {
      *(*sequence + i)  = hap_seq[i];
      *(*count + i) = hap_count[i];
    }
  *no_sequence = n_hapseq;
  regfree(second_regexCompiled);
  free(second_regexCompiled);
  free(pmatch);
  free(hap_seq);
  free(hap_count);
}

void iterator(gpointer key, gpointer value, gpointer user_data)  {
  printf(user_data, *(guint32*) key, *(guint32*)value);
}

void read_pop(char *filename, GIOChannel *infile, chrom_data *head, char pop_name)
{
  chrom_data* current;
  current = head;
  char *complete_file;
  if(g_io_channel_read_to_end(infile,&complete_file,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *source;
  source = complete_file;
  char * regexString = "([0-9]+)[[:blank:]]+(.+)";
  size_t maxGroups = 3;
                                                                                
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];

  if (regcomp(&regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    };

  /* printf("Reading haplotype information of population %c from file: %s ...\n", pop_name, filename); /\* reading file status printout *\/ */

  /* if(verbose) */
  /*   { */
  /*     printf("population | chrom index | #haplotypes\n"); */
  /*     printf("%.*s\n",38,"--------------------------------------------------------------------------------"); */
  /*   } */
  unsigned int chrom_no = 0;
  unsigned int *temporary_key;
  guint* temporary_count;
  while(current != NULL)
    {
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
	{
	  unsigned int g = 0;
	  unsigned int offset = 0;
	  char sourceCopy[strlen(source) + 1];
	  strcpy(sourceCopy, source);
	  sourceCopy[groupArray[g].rm_eo] = 0;
	  for (g = 0; g < maxGroups; g++)
	    {
	      if (groupArray[g].rm_so == (size_t)-1)
		break;  // No more groups

	      if(g == 0)
		{
		  offset = groupArray[g].rm_eo;
		}
	      if(g == 1)
		{
		  char *submatch_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
#ifdef DEBUG
		  printf("%s\n\n",submatch_1);
#endif
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  if(pop_name == 'A')
		    {		 
		      char *submatch_2;
		      size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_2 = (char*)malloc(matchlen_2+1);
		      assert(submatch_2 != NULL);
		      strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		      submatch_2[matchlen_2]='\0';
#ifdef DEBUG
		      printf("~~ %s ~~\n\n", submatch_2);
#endif
		      unsigned int *hap_sequence, no_hap_sequence;
		      unsigned int *hap_count;
		      hap_info_extract(submatch_2, &hap_sequence, &hap_count, &no_hap_sequence);
#ifdef DEBUG
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  printf("Sequence: %u\tCount: %u\n",hap_sequence[i], hap_count[i]);
			}
#endif
		      current->hashA = g_hash_table_new/* _full */(g_int_hash, g_int_equal/* , free, free */);
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  temporary_key = g_new(guint,1);
			  *temporary_key = hap_sequence[i];
			  temporary_count = g_new(guint,1);
			  *temporary_count = hap_count[i];
			  g_hash_table_insert(current->hashA, temporary_key, temporary_count);
			}
		      /* if(verbose) */
		      /* 	{ */
		      /* 	  printf("%10c | %11d | %11d\n", pop_name, current->chrom_id, g_hash_table_size(current->hashA)); */
		      /* 	} */
		      free(hap_sequence);
		      free(hap_count);
		      free(submatch_2);
		      
		    }
		  else
		    {
		      char *submatch_2;
		      size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_2 = (char*)malloc(matchlen_2+1);
		      assert(submatch_2 != NULL);
		      strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		      submatch_2[matchlen_2]='\0';
#ifdef DEBUG
		      printf("~~ %s ~~\n\n", submatch_2);
#endif
		      unsigned int *hap_sequence, no_hap_sequence;
		      unsigned int *hap_count;
		      hap_info_extract(submatch_2, &hap_sequence, &hap_count, &no_hap_sequence);
#ifdef DEBUG
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  printf("Sequence: %u\tCount: %u\n",hap_sequence[i], hap_count[i]);
			}
#endif
		      current->hashB = g_hash_table_new/* _full */(g_int_hash, g_int_equal/* , free, free */);
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  temporary_key = g_new(guint,1);
			  *temporary_key = hap_sequence[i];
			  temporary_count = g_new(guint,1);
			  *temporary_count = hap_count[i];
			  g_hash_table_insert(current->hashB, temporary_key, temporary_count);
			}
		      
		      /* printf("chromosome %d has %d haplotypes\n", current->chrom_id, g_hash_table_size(current->hashB)); */
		      /* if(verbose) */
		      /* 	{ */
		      /* 	  printf("%10c | %11d | %11d\n", pop_name, current->chrom_id, g_hash_table_size(current->hashB)); */
		      /* 	} */
		      free(hap_sequence);
		      free(hap_count);
		      free(submatch_2);

		    }

		}
	    }
	  source = source + offset;
	  chrom_no++;
	  /* printf("\n\n"); */
	}
      current = current->next;
    }
  /* printf("Finished reading the file containing haplotype information for population %c successfully!\n\n",pop_name); */
  regfree(&regexCompiled);
  g_free(complete_file);
}


void update_haplotype_indv(char *SNP, unsigned int *hap1, unsigned int *hap2, unsigned int position)
{
  char *second_regexString = "(0|1)\\|(0|1)";
  char *second_string;
  second_string = SNP;
/* #ifdef DEBUG */
/*   printf("%s\n\n",second_string); */
/* #endif */
  size_t max_subgroup = 3;
  regex_t second_regexCompiled;
  regmatch_t pmatch[max_subgroup];
  unsigned int m = 0;
  char *cursor;
  
  if(regcomp(&second_regexCompiled, second_regexString, REG_EXTENDED))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    }
  cursor = second_string;
  /* unsigned int flag; */


  while(regexec(&second_regexCompiled, cursor, max_subgroup, pmatch, 0) == 0)
    {

      /* flag = regexec(&second_regexCompiled, cursor, max_subgroup, pmatch, 0); */
      /* printf("~~ Flag:%u ~~\n",flag); */

      unsigned int g1 = 0;
      unsigned int offset_1 = 0;
      for (g1 = 0; g1 < max_subgroup; g1++)
	{
	  if (pmatch[g1].rm_so == (size_t)-1)
	    break;  // No more groups

	  if (g1 == 0)
	    offset_1 = pmatch[g1].rm_eo;

	  char cursorCopy[strlen(cursor) + 1];
	  strcpy(cursorCopy, cursor);
	  cursorCopy[pmatch[g1].rm_eo] = 0;
/* #ifdef DEBUG */
/* 	  printf("Match %u, Group %u: [%2u-%2u]: %s\n", */
/* 		 m, g1, pmatch[g1].rm_so, pmatch[g1].rm_eo,		  */
/* 		 cursorCopy + pmatch[g1].rm_so); */
/* #endif */
	  if(g1 == 1)
	  	{
	  	  unsigned int val_1;
	  	  val_1 = (unsigned int) atoi(cursorCopy + pmatch[g1].rm_so);
		  /* printf("%u\n",val_1); */
	  	  hap1[m] = (hap1[m] |(val_1 << position));
	  	}
	  if(g1 == 2)
	  	{
	  	  unsigned int val_2;
	  	  val_2 = (unsigned int) atoi(cursorCopy + pmatch[g1].rm_so);
		  /* printf("%u\n",val_2); */
	  	  hap2[m] = (hap2[m] |(val_2 << position));
	  	}

	}
      cursor += offset_1;
      m++;
    }

  regfree(&second_regexCompiled);

}


void update_linked_list_indv_data(char *individual_data, chrom_data *head, unsigned int no_indv)
{
  chrom_data *current;
  current = head;
  /* printf("~~~~\n\n"); */
  /* printf("%s\n", individual_data); */
  char *source;
  source = individual_data;
  char * regexString = "([0-9]+):([0-9]+)[[:blank:]]+(.+)";
  size_t maxGroups = 4;
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];

  if (regcomp(&regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    };
  
  while(current != NULL)
    {
      unsigned int *hap_1, *hap_2;
      hap_1 = calloc(no_indv, sizeof(unsigned int));
      hap_2 = calloc(no_indv, sizeof(unsigned int));
      assert(hap_1 != NULL);
      assert(hap_2 != NULL);
      for(int i = 0; i < current->n_loci; i++)
	{
	  if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
	    {
	      unsigned int g = 0;
	      unsigned int offset = 0;
	      char sourceCopy[strlen(source) + 1];
	      strcpy(sourceCopy, source);
	      sourceCopy[groupArray[g].rm_eo] = 0;

	      for (g = 0; g < maxGroups; g++)
		{
		  if (groupArray[g].rm_so == (size_t)-1)
		    break;  // No more groups

		  if(g == 0)
		    {
		      offset = groupArray[g].rm_eo;
		    }
		  if(g == 1)
		    {
		      char *submatch_1;
		      size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_1 = (char*)malloc(matchlen_1+1);
		      assert(submatch_1 != NULL);
		      strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		      submatch_1[matchlen_1]='\0';
		      /* printf("%s\n\n",submatch_1); */
		      free(submatch_1);
		    }
		  if(g == 2)
		    {
		      char *submatch_2;
		      size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_2 = (char*)malloc(matchlen_2+1);
		      assert(submatch_2 != NULL);
		      strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		      submatch_2[matchlen_2]='\0';
		      /* printf("~~ %s ~~\n\n", submatch_2); */
		      free(submatch_2);
		    }
		  if(g == 3)
		    {
		      char *submatch_3;
		      size_t matchlen_3 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_3 = (char*)malloc(matchlen_3+1);
		      assert(submatch_3 != NULL);
		      strncpy(submatch_3, sourceCopy + groupArray[g].rm_so, matchlen_3+1);
		      submatch_3[matchlen_3]='\0';
		      /* printf("~~ %s ~~\n\n", submatch_3); */
		      update_haplotype_indv(submatch_3, hap_1, hap_2, current->n_loci-i-1);
		      free(submatch_3);
		    }
		}
	      source = source + offset;
	    }
	}
      current->haplotype_1 = malloc(no_indv * sizeof(unsigned int));
      current->haplotype_2 = malloc(no_indv * sizeof(unsigned int));
      assert(current->haplotype_1 != NULL);
      assert(current->haplotype_2 != NULL);
      for(int j = 0; j < no_indv; j++)
	{
	  current->haplotype_1[j] = hap_1[j];
	  current->haplotype_2[j] = hap_2[j];
#ifdef DEBUG
	  printf("%u\t%u\n",hap_1[j],hap_2[j]);
#endif
	}
      /* printf("\n\n"); */
      free(hap_1);
      free(hap_2);
      current = current->next;
    }
  /* printf("Number of individuals in the data file: %u\n", no_indv); */
  /* printf("Finished reading the data file successfully!\n\n"); */
  regfree(&regexCompiled);

}

void read_multiple_indv(char *filename, GIOChannel *infile, chrom_data *head, unsigned int *n_indv/* , char *indv_names */)
{
  char *header;
  if(g_io_channel_read_line(infile,&header,NULL,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }

  char *regexString = "([[:alnum:]:_]+)";
  size_t max_subgroup = 1;
  regex_t regexCompiled;
  regmatch_t pmatch[max_subgroup];
  unsigned int m = 0;
  char *cursor;

  /* char *individual_ID[MAX_INDIVIDUALS]; */
  
  if(regcomp(&regexCompiled, regexString, REG_EXTENDED))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    }
  cursor = header;

  /* printf("Reading data file: %s ...\n", filename); */
  /* unsigned int flag; */
  while(regexec(&regexCompiled, cursor, max_subgroup, pmatch, 0) == 0)
    {

      /* flag = regexec(&regexCompiled, cursor, max_subgroup, pmatch, 0); */
      /* printf("~~ Flag:%u ~~\n",flag); */

      unsigned int g1 = 0;
      unsigned int offset = 0;
      for (g1 = 0; g1 < max_subgroup; g1++)
        {
          if (pmatch[g1].rm_so == (size_t)-1)
            break;  // No more groups

          if (g1 == 0)
            offset = pmatch[g1].rm_eo;

          char cursorCopy[strlen(cursor) + 1];
          strcpy(cursorCopy, cursor);
          cursorCopy[pmatch[g1].rm_eo] = 0;
          /* printf("Match %u, Group %u: [%2u-%2u]: %s\n", */
          /*        m, g1, pmatch[g1].rm_so, pmatch[g1].rm_eo, */
          /*        cursorCopy + pmatch[g1].rm_so); */
	  /* if(m != 0){ */
	  /*   strcpy(individual_ID[m],cursorCopy + pmatch[g1].rm_so); */
	  /* } */
        }
      cursor += offset;
      m++;
    }

  /* printf("Count regexec: %u\n",m); */


  
  
  char *complete_file_without_header;
  if(g_io_channel_read_to_end(infile,&complete_file_without_header,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *data;
  data = complete_file_without_header;
  /* printf("\n\n%s\n\n",data); */
  *n_indv = (m-1);

  update_linked_list_indv_data(data,head,m-1);

  regfree(&regexCompiled);
  g_free(complete_file_without_header);
  g_free(header);

}


void add_counts(gpointer key, gpointer value, gpointer userdata)
{
  unsigned int *u;
  u = userdata;

  unsigned int v = (unsigned int)*(guint*)value;
  *u += v;
}

GHashTable* update_hap(GHashTable *table, unsigned int *hap, unsigned int add)
{
  gpointer a;
  a = g_hash_table_lookup(table, hap);

  if(add == 1)
    {
      *(guint*) a = *(guint*) a + 1;
      g_hash_table_insert(table, hap, a);
    }
  else
    {
      *(guint*) a = *(guint*) a - 1;
      g_hash_table_insert(table, hap, a);

    }
  
  return(table);

}

double prob_change(double interval_length, double recom_rate)
{
  double lambda, value;
  
  lambda = (interval_length * recom_rate);
  value = 1 - (exp(-lambda) * cosh(lambda));
  return(value);
}

void all_marker_prob_change(double **inter_mark_prob_change,double **inter_marker_length, double recom_rate, unsigned int n_loci)
{
  for(int i = 0; i < n_loci; i++)
    {
      *(*inter_mark_prob_change + i) = prob_change(*(*inter_marker_length + i),recom_rate);
    }
}


unsigned int returnBit(unsigned int n, unsigned int k)
{
  unsigned int temp;
  /* temp = n; */
  temp = (n & (1 << k));
  if(temp == 0)
    {
      return(0);
    }
  else
    {
      return(1);
    }
}

double Q_component(unsigned int start_state, unsigned int ancestry_state,double **inter_marker_prob_change, unsigned int n_loci, double pi)
{
  unsigned int previous_state, present_state;
  previous_state = start_state;
  double log_product;
  log_product = log(pi);
  for(int i = 0; i < n_loci; i++)
    {
      present_state = returnBit(ancestry_state, n_loci - i -1);
      if(present_state != previous_state)
	{
	  log_product = log_product + log(*(*inter_marker_prob_change + i));
	}
      else
	{
	  log_product = log_product + log(1 - *(*inter_marker_prob_change + i));
	}
      previous_state = present_state;
    }
  double Q;
  Q = log_product;
  return(Q);
}

double Q(unsigned int ancestry_state, double **inter_marker_prob_change, unsigned int n_loci, double pi_0, double pi_1)
{
  double log_Q_0, log_Q_1;
  log_Q_0 = Q_component(0, ancestry_state, inter_marker_prob_change, n_loci, pi_0);
  log_Q_1 = Q_component(1, ancestry_state, inter_marker_prob_change, n_loci, pi_1);
  double scale;
  scale = fabs(fmax(log_Q_0,log_Q_1)) - 0.5;
  
  double log_Q;
  log_Q = log(exp(log_Q_0 + scale) + exp(log_Q_1 + scale)) - scale;
  return(log_Q);
  
}




unsigned int get_sub_hap(unsigned int ancestry_state, unsigned int hap, unsigned int n_loci, unsigned int population)
{
  unsigned int result, count, term;
  result = 0;
  count = 0;

  if (population == 0)
    {
      for(int i = 0; i < n_loci; i++)
	{                        
	  if (returnBit(ancestry_state,i) == 0)
	    {                    
	      term = (returnBit(hap,i) << count);
	      /* printf("%u\n", term); */
	      result = (result | term);
	      count++;
	    }
	}
    }
  else
    {
      for(int i = 0; i < n_loci; i++)
	{                        
	  if (returnBit(ancestry_state,i) == 1)
	    {                    
	      term = (returnBit(hap,i) << count);
	      /* printf("%u\n", term); */
	      result = (result | term);
	      count++;
	    }
	}

    }

  return(result);
}


double p_recom(unsigned int ancestry_state, unsigned int hap, unsigned int n_loci, GHashTable *table_A, GHashTable *table_B)
{
  
  GHashTable* hashA_copy = g_hash_table_new_full(g_int_hash, g_int_equal, NULL, NULL);
  GHashTable* hashB_copy = g_hash_table_new_full(g_int_hash, g_int_equal, NULL, NULL);

  unsigned int table_size_A, table_size_B;
  table_size_A = g_hash_table_size(table_A);
  table_size_B = g_hash_table_size(table_B);

  gpointer *keys_A, *keys_B;
  keys_A = g_hash_table_get_keys_as_array(table_A, &table_size_A);
  keys_B = g_hash_table_get_keys_as_array(table_B, &table_size_B);
  
  unsigned int *copy_key_A, *copy_value_A;
  copy_key_A = g_new(guint, table_size_A);
  copy_value_A = g_new(guint, table_size_A);

  copy_key_A[0] = get_sub_hap(ancestry_state, *(guint*) keys_A[0], n_loci, 0);
  copy_value_A[0] = *(guint*) g_hash_table_lookup(table_A, keys_A[0]);
  g_hash_table_insert(hashA_copy, &copy_key_A[0], &copy_value_A[0]);

  for(int i = 1; i < table_size_A; i++)
    {
      copy_key_A[i] = get_sub_hap(ancestry_state, *(guint*) keys_A[i], n_loci, 0);
      /* printf("*** %u \t%u ***\n", *(guint*) keys_A[i], copy_key_A[i]); */
      /* printf("CHECKPOST--------------------------->\n"); */

      if (g_hash_table_lookup(hashA_copy, &copy_key_A[i]) == NULL)
      	{
      	  copy_value_A[i] = *(guint*) g_hash_table_lookup(table_A, keys_A[i]);
      	  g_hash_table_insert(hashA_copy, &copy_key_A[i], &copy_value_A[i]);
      	}
      else
      	{
      	  copy_value_A[i] = *(guint*) g_hash_table_lookup(table_A, keys_A[i]) + *(guint*) g_hash_table_lookup(hashA_copy, &copy_key_A[i]);
      	  g_hash_table_insert(hashA_copy, &copy_key_A[i], &copy_value_A[i]);
      	}

    }

  /* g_hash_table_foreach(hashA_copy, (GHFunc)iterator, "Haplotype Sequence: %u has Counts %u \n"); */
  /* printf("*********************************************************************\n\n"); */
  
  unsigned int *copy_key_B, *copy_value_B;
  copy_key_B = g_new(guint, table_size_B);
  copy_value_B = g_new(guint, table_size_B);

  copy_key_B[0] = get_sub_hap(ancestry_state, *(guint*) keys_B[0], n_loci, 1);
  copy_value_B[0] = *(guint*) g_hash_table_lookup(table_B, keys_B[0]);
  g_hash_table_insert(hashB_copy, &copy_key_B[0], &copy_value_B[0]);


  for(int i = 1; i < table_size_B; i++)
    {
      copy_key_B[i] = get_sub_hap(ancestry_state, *(guint*) keys_B[i], n_loci, 1);
      /* printf("* %u -- %u *\n", *(guint*) keys_B[i], copy_key_B[i]); */
      
      if (g_hash_table_lookup(hashB_copy, &copy_key_B[i]) == NULL)
      	{
      	  copy_value_B[i] = *(guint*) g_hash_table_lookup(table_B, keys_B[i]);
      	  g_hash_table_insert(hashB_copy, &copy_key_B[i], &copy_value_B[i]);
      	}
      else
      	{
      	  copy_value_B[i] = *(guint*) g_hash_table_lookup(table_B, keys_B[i]) + *(guint*) g_hash_table_lookup(hashB_copy, &copy_key_B[i]);
      	  g_hash_table_insert(hashB_copy, &copy_key_B[i], &copy_value_B[i]);
      	}

    }

  /* g_hash_table_foreach(hashB_copy, (GHFunc)iterator, "Haplotype Sequence: %u has Counts %u \n"); */

  unsigned int *total_countA, *total_countB;
  total_countA = g_new(guint,1);
  total_countB = g_new(guint,1);
  *total_countA = 0;
  *total_countB = 0;
  g_hash_table_foreach(hashA_copy, add_counts, (gpointer) total_countA);
  g_hash_table_foreach(hashB_copy, add_counts, (gpointer) total_countB);

  unsigned int sub_table_size_A, sub_table_size_B;
  sub_table_size_A = g_hash_table_size(hashA_copy);
  sub_table_size_B = g_hash_table_size(hashB_copy);

  
  unsigned int sub_hap_A, sub_hap_B;

  
  double value;
  unsigned int all_B_anc_state;
  all_B_anc_state = pow(2,n_loci) - 1;
  if (ancestry_state == 0 || ancestry_state == all_B_anc_state)
    {
      if (ancestry_state == 0)
  	{
  	  /* update_hap(hashA_copy, &hap, add_hap); */
  	  /* value = p_prob(hashA_copy); */
	  unsigned int n_hapA;
	  if(g_hash_table_lookup(hashA_copy, &hap) == NULL)
	    {
	      n_hapA = 0;
	    }
	  else
	    {
	      n_hapA = *(guint*) g_hash_table_lookup(hashA_copy, &hap);
	    }
	  value = lgamma(1+ *total_countA) - lgamma(2+ *total_countA)
	    + lgamma(1+(1/(double)sub_table_size_A)+n_hapA) - lgamma((1/(double)sub_table_size_A)+n_hapA);
  	}
      else
  	{
  	  /* update_hap(hashB_copy, &hap, add_hap); */
  	  /* value = p_prob(hashB_copy); */
	  unsigned int n_hapB;
	  if(g_hash_table_lookup(hashB_copy, &hap) == NULL)
	    {
	      n_hapB = 0;
	    }
	  else
	    {
	      n_hapB = *(guint*) g_hash_table_lookup(hashB_copy, &hap);
	    }
	  value = lgamma(1+ *total_countB) - lgamma(2+ *total_countB)
	    + lgamma(1+(1/(double)sub_table_size_B)+n_hapB) - lgamma((1/(double)sub_table_size_B)+n_hapB);
  	}
    }
  else
    {
      sub_hap_A = get_sub_hap(ancestry_state, hap, n_loci, 0);
      sub_hap_B = get_sub_hap(ancestry_state, hap, n_loci, 1);
      /* update_hap(hashA_copy, &sub_hap_A, add_hap); */
      /* update_hap(hashB_copy, &sub_hap_B, add_hap); */
      /* value = p_prob(hashA_copy) + p_prob(hashB_copy); */

      unsigned int n_hapA, n_hapB;
      if(g_hash_table_lookup(hashA_copy, &sub_hap_A) == NULL)
	{
	  n_hapA = 0;
	}
      else
	{
	  n_hapA = *(guint*) g_hash_table_lookup(hashA_copy, &sub_hap_A);
	}

      if(g_hash_table_lookup(hashB_copy, &sub_hap_B) == NULL)
	{
	  n_hapB = 0;
	}
      else
	{
	  n_hapB = *(guint*) g_hash_table_lookup(hashB_copy, &sub_hap_B);
	}
      
      value = lgamma(1+ *total_countA) - lgamma(2+ *total_countA)
	+ lgamma(1+(1/(double)sub_table_size_A)+n_hapA) - lgamma((1/(double)sub_table_size_A)+n_hapA)
	+ lgamma(1+ *total_countB) - lgamma(2+ *total_countB)
	+ lgamma(1+(1/(double)sub_table_size_B)+n_hapB) - lgamma((1/(double)sub_table_size_B)+n_hapB);

      
    }

  g_free(total_countA);
  g_free(total_countB);
  
  g_free(copy_key_A);
  g_free(copy_value_A);
  g_hash_table_destroy(hashA_copy);
  g_free(keys_A);

  g_free(copy_key_B);
  g_free(copy_value_B);
  g_hash_table_destroy(hashB_copy);
  g_free(keys_B);

  return(value);
}

double recom_prob(unsigned int hap, chrom_data *data_node)
{
  double *loci_position;
  loci_position = malloc(data_node->n_loci * sizeof(double));
  assert(loci_position != NULL);
  for(int i = 0; i < data_node->n_loci; i++)
    {
      loci_position[i] = (double) (data_node->markers[i] / (data_node->chrom_length * pow(10,6)));
    }

  double *interval_length;
  interval_length = malloc(data_node->n_loci * sizeof(double));
  assert(interval_length != NULL);
  for(int i = 0; i < data_node->n_loci; i++)
    {
      if(i == 0)
  	{
  	  interval_length[i] = loci_position[i];
  	}
      else
	{
	  interval_length[i] = loci_position[i] - loci_position[i-1];
	}
    }

  double r; /* r is the expectated number of recombination for a given length of chromosome */
  r = (data_node->chrom_recom_rate * data_node->chrom_length / 100.0);
  
  /* printf("The probability of state change for each inter-marker interval is given by\n"); */
  double *P;
  P = malloc(data_node->n_loci * sizeof(double));
  assert(P != NULL);
  all_marker_prob_change(&P,&interval_length,r,data_node->n_loci);
/* #ifdef DEBUG */
/*   for(int i = 0; i < data_node->n_loci; i++) */
/*     { */
/*       printf("P[%d] = %lf\n",i,P[i]); */
/*     } */
/* #endif */
  
  /* printf("\n\n"); */
  /* printf("The value of Q for each ancestry state is as follows:\n"); */
  unsigned int index, mid_index;
  index = pow(2,data_node->n_loci);
  mid_index = pow(2,data_node->n_loci - 1);
  double *Q_ancestry;
  Q_ancestry = malloc(index * sizeof(double));
  assert(Q_ancestry != NULL);
  for(unsigned int i = 0; i < mid_index; i++)
    {
      Q_ancestry[i] = Q(i, &P, data_node->n_loci, 0.5, 0.5);
      Q_ancestry[index-i-1] = Q_ancestry[i];
    }

/* #ifdef DEBUG */
/*   for(unsigned int i = 0; i < index; i++) */
/*     { */
/*       printf("Q[%u] = %lf\n",i,Q_ancestry[i]); */
/*     } */
/* #endif */


  /* chrom_data* data_node_copy; */
  /* data_node_copy = data_node; */

  double scale, maximum, temporary;
  double *log_p, *log_R;
  log_p = malloc(index * sizeof(double));
  log_R = malloc(index * sizeof(double));
  assert(log_p != NULL);
  assert(log_R != NULL);
  double final_log_prob;
  double sum = 0;
  
  log_p[0] = p_recom(0,hap, data_node->n_loci, data_node->hashA, data_node->hashB);
  /* printf("### %lf\n", log_p[0]); */
  log_R[0] = log_p[0] + Q_ancestry[0];
  maximum = log_R[0];
  temporary = log_R[0];
  for(unsigned int i = 1; i < index; i++)
    {
      log_p[i] = p_recom(i,hap,data_node->n_loci,data_node->hashA,data_node->hashB);

      /* printf("### %lf\n", log_p[i]); */
      log_R[i] = log_p[i] + Q_ancestry[i];

      temporary = log_R[i];
      if(maximum < temporary)
	{
	  maximum = temporary;
	}
    }
  /* printf("%lf\n",maximum); */
  scale = fabs(maximum)-0.5;
  for(unsigned int i = 0; i < index; i++)
    {
      sum = sum + exp(log_R[i] + scale) ;
    }
  final_log_prob = log(sum) - scale;
  /* printf("\n\nLog(base e) of the recombinant haplotype: %lf\n\n", final_log_prob); */
  free(log_p);
  free(log_R);


  free(interval_length);
  free(P);
  free(Q_ancestry);
  free(loci_position);

  return(final_log_prob);
}



double p_model(chrom_data* node, unsigned int hap1, unsigned int hap2, unsigned int model)
  /* double p_model(GHashTable *table_A, GHashTable *table_B, unsigned int hap1, unsigned int hap2, unsigned int model) */
{
  double value;
  GHashTable* hashA_copy = g_hash_table_new_full(g_int_hash, g_int_equal, NULL, NULL);
  GHashTable* hashB_copy = g_hash_table_new_full(g_int_hash, g_int_equal, NULL, NULL);

  unsigned int table_size_A, table_size_B;
  table_size_A = g_hash_table_size(node->hashA);
  table_size_B = g_hash_table_size(node->hashB);

  gpointer *keys_A, *keys_B;
  keys_A = g_hash_table_get_keys_as_array(node->hashA, &table_size_A);
  keys_B = g_hash_table_get_keys_as_array(node->hashB, &table_size_B);
  

  unsigned int *copy_key_A, *copy_value_A;
  copy_key_A = g_new(guint, table_size_A);
  copy_value_A = g_new(guint, table_size_A);

  unsigned int total_hap_A;
  total_hap_A = 0;
  for(int i = 0; i < table_size_A; i++)
    {      
      copy_key_A[i] = *(guint*) keys_A[i];
      copy_value_A[i] = *(guint*) g_hash_table_lookup(node->hashA, keys_A[i]);
      total_hap_A = total_hap_A + copy_value_A[i];
      g_hash_table_insert(hashA_copy, &copy_key_A[i], &copy_value_A[i]);
    }
  
  unsigned int *copy_key_B, *copy_value_B;
  copy_key_B = g_new(guint, table_size_B);
  copy_value_B = g_new(guint, table_size_B);

  unsigned int total_hap_B;
  total_hap_B = 0;

  for(int i = 0; i < table_size_B; i++)
    {      
      copy_key_B[i] = *(guint*) keys_B[i];
      copy_value_B[i] = *(guint*) g_hash_table_lookup(node->hashB, keys_B[i]);
      total_hap_B = total_hap_B + copy_value_B[i];
      g_hash_table_insert(hashB_copy, &copy_key_B[i], &copy_value_B[i]);
    }

  /* printf("Total #hap in A: %u\n", total_hap_A); */
  /* printf("Total #hap in B: %u\n", total_hap_B); */
  
  /* printf("#distinct hap in A: %u\n", table_size_A); */
  /* printf("#distinct hap in B: %u\n", table_size_B); */

  /* unsigned int *total_count; */
  /* total_count = g_new(guint,1); */
  /* *total_count = 0; */
  /* g_hash_table_foreach(hashB_copy, add_counts, (gpointer) total_count); */

  /* printf("~~~~~~~~~~ Total #hap in A: %u\n", *total_count); */
  /* g_free(total_count); */
  
  if (model == 0)
    {
      /* update_hap(hashB_copy, &hap1, add_hap); */
      /* update_hap(hashB_copy, &hap2, add_hap); */
      /* value = p_prob(hashA_copy) + p_prob(hashB_copy); */

      if(hap1 == hap2)
	{
	  unsigned int n_hap1;
	  if(g_hash_table_lookup(hashB_copy, &hap1) == NULL)
	    {
	      n_hap1 = 0;
	    }
	  else
	    {
	      n_hap1 = *(guint*) g_hash_table_lookup(hashB_copy, &hap1);
	    }
	  value = lgamma(1+total_hap_B) - lgamma(3+total_hap_B)
	    + lgamma(2+(1/(double)table_size_B)+n_hap1) - lgamma((1/(double)table_size_B)+n_hap1);
	}
      else
	{
	  unsigned int n_hap1, n_hap2;
	  if(g_hash_table_lookup(hashB_copy, &hap1) == NULL)
	    {
	      n_hap1 = 0;
	    }
	  else
	    {
	      n_hap1 = *(guint*) g_hash_table_lookup(hashB_copy, &hap1);
	    }

	  if(g_hash_table_lookup(hashB_copy, &hap2) == NULL)
	    {
	      n_hap2 = 0;
	    }
	  else
	    {
	      n_hap2 = *(guint*) g_hash_table_lookup(hashB_copy, &hap2);
	    }
	  value = lgamma(1+total_hap_B) - lgamma(3+total_hap_B) + lgamma(3) - 2*lgamma(2)
	    + lgamma(1+(1/(double)table_size_B)+n_hap1) - lgamma((1/(double)table_size_B)+n_hap1)
	    + lgamma(1+(1/(double)table_size_B)+n_hap2) - lgamma((1/(double)table_size_B)+n_hap2);
	}
      
      /* printf("Model (a): %lf\n\n", value); */
      
    }

  else if (model == 1)
    {
      if (hap1 == hap2)
	{
	  unsigned int n_hap1_A;
	  if(g_hash_table_lookup(hashA_copy, &hap1) == NULL)
	    {
	      n_hap1_A = 0;
	    }
	  else
	    {
	      n_hap1_A = *(guint*) g_hash_table_lookup(hashA_copy, &hap1);
	    }
	  value = lgamma(1+total_hap_A) - lgamma(2+total_hap_A)
	    + lgamma(1+(1/(double)table_size_A)+n_hap1_A) - lgamma((1/(double)table_size_A)+n_hap1_A)
	    + recom_prob(hap2, node);
	}
      else
	{
	  double term1_b, term2_b;
	  /* update_hap(hashA_copy, &hap1, add_hap); */

	  unsigned int n_hap1_A, n_hap2_A;
	  if(g_hash_table_lookup(hashA_copy, &hap1) == NULL)
	    {
	      n_hap1_A = 0;
	    }
	  else
	    {
	      n_hap1_A = *(guint*) g_hash_table_lookup(hashA_copy, &hap1);
	    }
	  term1_b = lgamma(1+total_hap_A) - lgamma(2+total_hap_A)
	    + lgamma(1+(1/(double)table_size_A)+n_hap1_A) - lgamma((1/(double)table_size_A)+n_hap1_A)
	    + recom_prob(hap2, node);
      
	  /* term1_b = p_prob(hashA_copy) + recom_prob(hap2, node);  */
	  /* update_hap(hashA_copy, &hap1, remove_hap); */
	  /* update_hap(hashA_copy, &hap2, add_hap); */

	  if(g_hash_table_lookup(hashA_copy, &hap2) == NULL)
	    {
	      n_hap2_A = 0;
	    }
	  else
	    {
	      n_hap2_A = *(guint*) g_hash_table_lookup(hashA_copy, &hap2);
	    }
	  term2_b = lgamma(1+total_hap_A) - lgamma(2+total_hap_A)
	    + lgamma(1+(1/(double)table_size_A)+n_hap2_A) - lgamma((1/(double)table_size_A)+n_hap2_A)
	    + recom_prob(hap1, node);
	  /* term2_b = p_prob(hashA_copy) + recom_prob(hap1, node); */


	  double scale_b;       /* scaling factors to avoid underflows */
	  scale_b = fabs(fmax(term1_b, term2_b)) - 0.5;
	  value = log(exp(term1_b + scale_b) + exp(term2_b + scale_b)) - scale_b;
	}

      /* printf("Model (b): %lf\n\n", value); */
    }
  else if (model == 2)
    {
      if (hap1 == hap2)
	{
	  unsigned int n_hap1_A, n_hap2_B;
	  if(g_hash_table_lookup(hashA_copy, &hap1) == NULL)
	    {
	      n_hap1_A = 0;
	    }
	  else
	    {
	      n_hap1_A = *(guint*) g_hash_table_lookup(hashA_copy, &hap1);
	    }

	  if(g_hash_table_lookup(hashB_copy, &hap2) == NULL)
	    {
	      n_hap2_B = 0;
	    }
	  else
	    {
	      n_hap2_B = *(guint*) g_hash_table_lookup(hashB_copy, &hap2);
	    }
	  value = lgamma(1+total_hap_A) - lgamma(2+total_hap_A)
	    + lgamma(1+(1/(double)table_size_A)+n_hap1_A) - lgamma((1/(double)table_size_A)+n_hap1_A)
	    + lgamma(1+total_hap_B) - lgamma(2+total_hap_B)
	    + lgamma(1+(1/(double)table_size_B)+n_hap2_B) - lgamma((1/(double)table_size_B)+n_hap2_B);
	}
      else
	{
	  /* update_hap(hashA_copy, &hap1, add_hap); */
	  /* update_hap(hashB_copy, &hap2, add_hap); */
	  double term1, term2;
	  /* term1 = p_prob(hashA_copy) + p_prob(hashB_copy); */
	  /* update_hap(hashA_copy, &hap1, remove_hap); */
	  /* update_hap(hashB_copy, &hap2, remove_hap); */

	  /* update_hap(hashA_copy, &hap2, add_hap); */
	  /* update_hap(hashB_copy, &hap1, add_hap); */
	  /* term2 = p_prob(hashA_copy) + p_prob(hashB_copy); */

	  unsigned int n_hap1_A, n_hap2_A, n_hap1_B, n_hap2_B;
	  if(g_hash_table_lookup(hashA_copy, &hap1) == NULL)
	    {
	      n_hap1_A = 0;
	    }
	  else
	    {
	      n_hap1_A = *(guint*) g_hash_table_lookup(hashA_copy, &hap1);
	    }
	  
	  if(g_hash_table_lookup(hashA_copy, &hap2) == NULL)
	    {
	      n_hap2_A = 0;
	    }
	  else
	    {
	      n_hap2_A = *(guint*) g_hash_table_lookup(hashA_copy, &hap2);
	    }

	  if(g_hash_table_lookup(hashB_copy, &hap1) == NULL)
	    {
	      n_hap1_B = 0;
	    }
	  else
	    {
	      n_hap1_B = *(guint*) g_hash_table_lookup(hashB_copy, &hap1);
	    }

	  if(g_hash_table_lookup(hashB_copy, &hap2) == NULL)
	    {
	      n_hap2_B = 0;
	    }
	  else
	    {
	      n_hap2_B = *(guint*) g_hash_table_lookup(hashB_copy, &hap2);
	    }

	  term1 = lgamma(1+total_hap_A) - lgamma(2+total_hap_A)
	    + lgamma(1+(1/(double)table_size_A)+n_hap1_A) - lgamma((1/(double)table_size_A)+n_hap1_A)
	    + lgamma(1+total_hap_B) - lgamma(2+total_hap_B)
	    + lgamma(1+(1/(double)table_size_B)+n_hap2_B) - lgamma((1/(double)table_size_B)+n_hap2_B);

	  term2 = lgamma(1+total_hap_A) - lgamma(2+total_hap_A)
	    + lgamma(1+(1/(double)table_size_A)+n_hap2_A) - lgamma((1/(double)table_size_A)+n_hap2_A)
	    + lgamma(1+total_hap_B) - lgamma(2+total_hap_B)
	    + lgamma(1+(1/(double)table_size_B)+n_hap1_B) - lgamma((1/(double)table_size_B)+n_hap1_B);

      
	  double scale;
	  scale = fabs(fmax(term1, term2)) - 0.5;
	  value = log(exp(term1 + scale) + exp(term2 + scale)) - scale;
	}

      /* printf("Model (c): %lf\n\n", value); */
    }

  else if (model == 3)
    {
      /* update_hap(hashA_copy, &hap1, add_hap); */
      /* update_hap(hashA_copy, &hap2, add_hap); */
      /* value = p_prob(hashA_copy) + p_prob(hashB_copy); */

      if(hap1 == hap2)
	{
	  unsigned int n_hap1;
	  if(g_hash_table_lookup(hashA_copy, &hap1) == NULL)
	    {
	      n_hap1 = 0;
	    }
	  else
	    {
	      n_hap1 = *(guint*) g_hash_table_lookup(hashA_copy, &hap1);
	    }
	  value = lgamma(1+total_hap_A) - lgamma(3+total_hap_A)
	    + lgamma(2+(1/(double)table_size_A)+n_hap1) - lgamma((1/(double)table_size_A)+n_hap1);
	}
      else
	{
	  unsigned int n_hap1, n_hap2;
	  if(g_hash_table_lookup(hashA_copy, &hap1) == NULL)
	    {
	      n_hap1 = 0;
	    }
	  else
	    {
	      n_hap1 = *(guint*) g_hash_table_lookup(hashA_copy, &hap1);
	    }

	  if(g_hash_table_lookup(hashA_copy, &hap2) == NULL)
	    {
	      n_hap2 = 0;
	    }
	  else
	    {
	      n_hap2 = *(guint*) g_hash_table_lookup(hashA_copy, &hap2);
	    }

	  value = lgamma(1+total_hap_A) - lgamma(3+total_hap_A) + lgamma(3) - 2*lgamma(2)
	    + lgamma(1+(1/(double)table_size_A)+n_hap1) - lgamma((1/(double)table_size_A)+n_hap1)
	    + lgamma(1+(1/(double)table_size_A)+n_hap2) - lgamma((1/(double)table_size_A)+n_hap2);
	}

      /* printf("Model (d): %lf\n\n", value); */

    }
  
  else if (model == 4)
    {
      if(hap1 == hap2)
	{
	  unsigned int n_hap1_B;
	  if(g_hash_table_lookup(hashB_copy, &hap1) == NULL)
	    {
	      n_hap1_B = 0;
	    }
	  else
	    {
	      n_hap1_B = *(guint*) g_hash_table_lookup(hashB_copy, &hap1);
	    }
	  value = lgamma(1+total_hap_B) - lgamma(2+total_hap_B)
	    + lgamma(1+(1/(double)table_size_B)+n_hap1_B) - lgamma((1/(double)table_size_B)+n_hap1_B)
	    + recom_prob(hap2, node);
	}
      else
	{
	  double term1_e, term2_e;
	  /* update_hap(hashB_copy, &hap1, add_hap); */
	  unsigned int n_hap1_B, n_hap2_B;
	  if(g_hash_table_lookup(hashB_copy, &hap1) == NULL)
	    {
	      n_hap1_B = 0;
	    }
	  else
	    {
	      n_hap1_B = *(guint*) g_hash_table_lookup(hashB_copy, &hap1);
	    }

	  term1_e = lgamma(1+total_hap_B) - lgamma(2+total_hap_B)
	    + lgamma(1+(1/(double)table_size_B)+n_hap1_B) - lgamma((1/(double)table_size_B)+n_hap1_B)
	    + recom_prob(hap2, node);
      
	  /* term1_e = p_prob(hashB_copy) + recom_prob(hap2, node); */
	  /* update_hap(hashB_copy, &hap1, remove_hap); */

	  /* update_hap(hashB_copy, &hap2, add_hap); */
	  if(g_hash_table_lookup(hashB_copy, &hap2) == NULL)
	    {
	      n_hap2_B = 0;
	    }
	  else
	    {
	      n_hap2_B = *(guint*) g_hash_table_lookup(hashB_copy, &hap2);
	    }
	  term2_e = lgamma(1+total_hap_B) - lgamma(2+total_hap_B)
	    + lgamma(1+(1/(double)table_size_B)+n_hap2_B) - lgamma((1/(double)table_size_B)+n_hap2_B)
	    + recom_prob(hap1, node);
	  /* term2_e = p_prob(hashB_copy) + recom_prob(hap1, node); */

	  double scale_e;
	  scale_e = fabs(fmax(term1_e, term2_e)) - 0.5;
	  value = log(exp(term1_e + scale_e) + exp(term2_e + scale_e)) - scale_e;
	}

      /* printf("Model (e): %lf\n\n", value); */
    }

  else if (model == 5)
    {
      if(hap1 == hap2)
	{
	  value = recom_prob(hap1, node) + recom_prob(hap2, node);
	}
      else
	{
	  value = log(2) + recom_prob(hap1, node) + recom_prob(hap2, node);
	}
      /* printf("Model (f): %lf\n\n", value); */
    }

  
  g_free(copy_key_A);
  g_free(copy_value_A);
  g_hash_table_destroy(hashA_copy);
  g_free(keys_A);

  g_free(copy_key_B);
  g_free(copy_value_B);
  g_hash_table_destroy(hashB_copy);
  g_free(keys_B);
  
  return(value);
}


void largest(double array[], int n, GHashTable *hash_table,
	     char *best_model, char *next_best_model, double *bayes_factor)
{
  double largest, second_largest;
  int max_index, second_max_index;
  
  // Initialize largest and second largest element

  if(array[0] > array[1])
    {
      largest = array[0];
      second_largest = array[1];
      max_index = 0;
      second_max_index = 1;
    }
  else
    {
      largest = array[1];
      second_largest = array[0];
      max_index = 1;
      second_max_index = 0;
    }
  
  // Traverse array elements from third and
  // compare every element with current max and update second largest value
  for (int i = 2; i < n; i++)
    {
      if (array[i] > largest)
	{
	  second_largest = largest;
	  second_max_index = max_index;
	  largest = array[i];
	  max_index = i;
	}
      else
	{
	  if (array[i] > second_largest)
	    {
	      second_largest = array[i];
	      second_max_index = i;
	    }
	}
    }

  *best_model = *(gchar*)g_hash_table_lookup(hash_table,&max_index);
  *next_best_model = *(gchar*)g_hash_table_lookup(hash_table,&second_max_index);

  *bayes_factor = (largest/second_largest);
}



void free_a_hash_table_entry(gpointer key, gpointer value, gpointer user_data)
{
  g_free(key);
  g_free(value);
}

void free_all_key_value_entries (GHashTable *table)
{
    g_hash_table_foreach (table, free_a_hash_table_entry, NULL);
    g_hash_table_destroy (table);
}

void delete_chrom_data(chrom_data *head)
{
  chrom_data *temp;
  while(head != NULL)
    {
      temp = head;
      head = head->next;
      free(temp->markers);
      free_all_key_value_entries (temp->hashA);
      free_all_key_value_entries (temp->hashB);
      free(temp->haplotype_1);
      free(temp->haplotype_2);
      /* g_hash_table_destroy(temp->hashA); */
      /* g_hash_table_destroy(temp->hashB); */
      free(temp);
    }
}

int main(int argc, char **argv)
{
  
  int c;

  fillheader(version);
  show_header();
  
  int cflag = 0, Aflag = 0, Bflag = 0, iflag = 0, hflag = 0, Vflag = 0;
  char *filenames[4];
  char *output_filename = "out.txt";
  
  /* for(int i = 0; i < 4; i++) */
  /*   { */
  /*     filenames[i] = calloc(MAX_FILENAME, sizeof(char)); */
  /*     assert(filenames[i] != NULL); */
  /*   } */
  while((c = getopt(argc, argv, ":hVc:A:B:i:o:")) != -1)
    {
      switch(c)
	{
	case 'c':
	  /* strcpy(filenames[0],optarg); */
	  cflag = 1;
	  filenames[0] = optarg;
	  break;
	case 'A':
	  /* strcpy(filenames[1],optarg); */
	  Aflag = 1;
	  filenames[1] = optarg;
	  break;
	case 'B':
	  /* strcpy(filenames[2],optarg); */
	  Bflag = 1;
	  filenames[2] = optarg;
	  break;
	case 'i':
	  /* strcpy(filenames[3],optarg); */
	  iflag = 1;
	  filenames[3] = optarg;
	  break;
	case 'h':
	  hflag = 1;
	  break;
	case 'V':
	  Vflag = 1;
	  break;
	case 'o':
	  output_filename = optarg;
	  break;
	case ':':
	  if(optopt == 'c'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'A'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'B'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'i'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'o'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}

	  /* for(int i = 0; i < 4; i++) */
	  /*   { */
	  /*     free(filenames[i]); */
	  /*   } */
	  exit(1);
	case '?':
	  if (isprint (optopt))
	    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
	  else
	    fprintf(stderr,"Unknown option character `\\x%x'.\n",optopt);
	  /* for(int i = 0; i < 4; i++) */
	  /*   { */
	  /*     free(filenames[i]); */
	  /*   } */
	  exit(1);
	}
    }

  if(hflag)
    {
      program_help();
      exit(0);
    }

  if(optind < argc)
    {
      if(hflag)
	program_help();
      else
	print_msg();
      exit(1);
    }
  
  if(cflag == 0 || Aflag == 0 || Bflag == 0 || iflag == 0)
    {
      if(cflag == 0 && Aflag == 0 && Bflag == 0 && iflag == 0)
	{
	  if(hflag)
	    program_help();
	  else
	    print_msg();
	  exit(1);
	}
      if(cflag == 0)
	{
	  fprintf(stderr,"%s: missing -c option\n",argv[0]);
	}
      if(Aflag == 0)
	{
	  fprintf(stderr,"%s: missing -A option\n",argv[0]);
	}
      if(Bflag == 0)
	{
	  fprintf(stderr,"%s: missing -B option\n",argv[0]);
	}
      if(iflag == 0)
      	{
      	  fprintf(stderr,"%s: missing -i option\n",argv[0]);
      	}
      /* for(int i = 0; i < 4; i++) */
      /* 	{ */
      /* 	  free(filenames[i]); */
      /* 	} */
      exit(1);
    }

  /* get_main(cflag,Aflag,Bflag,iflag,filenames[0],filenames[1],filenames[2],filenames[3]); */
  
  if(cflag)
    {
      infile_chrom = g_io_channel_new_file(filenames[0],"r",&error);
      if(infile_chrom == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[0]);
	  exit(1);
	}
    }
  if(Aflag)
    {
      infile_A = g_io_channel_new_file(filenames[1],"r",&error);
      if(infile_A == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[1]);
	  exit(1);
	}

    }
  if(Bflag)
    {
      infile_B = g_io_channel_new_file(filenames[2],"r",&error);
      if(infile_B == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[2]);
	  exit(1);
	}
    }
  if(iflag)
    {
      infile_indv = g_io_channel_new_file(filenames[3],"r",&error);
      if(infile_indv == NULL)
  	{
  	  fprintf(stderr,"Could not open file %s\n", filenames[3]);
  	  exit(1);
  	}

    }


  unsigned int *chrom_index, *chrom_no_loci;
  unsigned int no_chrom;
  double *chrom_length, *chrom_recom_rate;
  double *loci;
  unsigned int total_markers;
  read_chrom(filenames[0],infile_chrom,&chrom_index,&chrom_no_loci,&chrom_length,&chrom_recom_rate,&no_chrom,&loci,&total_markers);

#ifdef DEBUG
  for(int i = 0; i < total_markers; i++)
    {
      printf("%lf\t",loci[i]);
    }
#endif
  chrom_data *head, *current;
  head = data_linked_list_creation(chrom_index, chrom_no_loci, chrom_length, chrom_recom_rate, no_chrom, loci);
  current = head;
  unsigned int n_chromosomes = 0;
  while(current != NULL)
    {
#ifdef DEBUG
      printf("chromosome number:%u\tchromosome length:%lf\trate_of_recom:%lf\tnumber_of_loci:%u\n",current->chrom_id,current->chrom_length, current->chrom_recom_rate, current->n_loci);
      for(int i = 0; i < current->n_loci; i++)
      	{
      	  printf("Marker[%d]:%lf\n",i+1,current->markers[i]);
      	}
#endif
      n_chromosomes++;
      current = current->next;
    }

  /* printf("\nNo. of chromosomes:%u\tTotal no.of markers:%u\n",n_chromosomes,total_markers); */
  char population_name[] = {'A','B','\0'};
  read_pop(filenames[1],infile_A,head,population_name[0]);

#ifdef DEBUG
  current = head;
  while(current != NULL)
    {
      g_hash_table_foreach(current->hashA, (GHFunc)iterator, "Haplotype Sequence: %u has Counts %u \n");
      printf("\n\n");
      current = current->next;
    }
#endif

  read_pop(filenames[2],infile_B,head,population_name[1]);

#ifdef DEBUG
  current = head;
  while(current != NULL)
    {
      g_hash_table_foreach(current->hashB, (GHFunc)iterator, "Haplotype Sequence: %u has Counts %u \n");
      printf("\n\n");
      current = current->next;
    }
#endif
  
  unsigned int individual_count;
  /* char *names_individuals; */
  /* char individual_name[] */
  read_multiple_indv(filenames[3], infile_indv, head, &individual_count/* , &names_individuals */);
  /* printf("%s\n",names_individuals); */
  
/* #ifdef DEBUG */
/*   current = head; */
/*   while(current != NULL) */
/*     { */
/*       printf("********* chromosome ********\n\n"); */
/*       for(int j = 0; j < individual_count; j++) */
/* 	{ */
/* 	  printf("%u\t%u\n",current->haplotype_1[j], current->haplotype_2[j]); */
/* 	  printf("%lf\n", p_model(current, current->haplotype_1[j], current->haplotype_2[j], 0)); */
/* 	  printf("%lf\n", p_model(current, current->haplotype_1[j], current->haplotype_2[j], 3)); */
/* 	  printf("%lf\n", p_model(current, current->haplotype_1[j], current->haplotype_2[j], 2)); */
/* 	  /\* printf("%lf\n", recom_prob(current->haplotype_2[j], current)); *\/ */
/* 	} */
/*       current = current->next; */
/*     } */
/* #endif */



  /* current = head; */
  /* while(current != NULL) */
  /*   { */
  /*     for(int j = 0; j < individual_count; j++) */
  /* 	{ */
  /* 	  printf("%u\t%u\n",current->haplotype_1[j], current->haplotype_2[j]); */
  /* 	} */
  /*     current = current->next; */
  /*   } */

  /* printf ("~~~ %lf ~~~~~\n", p_prob(head->hashA)); */
  /* printf ("~~~ %lf ~~~~~\n", p_prob(head->hashB)); */
  /* printf ("~~~ %lf ~~~~~\n", p_recom(1023,924,10,head->hashA,head->hashB)); */
  /* printf ("~~~ %lf ~~~~~\n", p_recom(783,924,10,head->hashA,head->hashB)); */


  /* -------------------------------------------------------------------------------------------------------------------------- */

  
  double *indv_model_prob[individual_count];
  for(int i = 0; i < individual_count; i++)
    {
      indv_model_prob[i] = (double*) calloc(NO_MODELS, sizeof(double));
      assert(indv_model_prob[i] != NULL);
    }

  double *final_prob[individual_count];
  for(int i = 0; i < individual_count; i++)
    {
      final_prob[i] = (double*) calloc(NO_MODELS, sizeof(double));
      assert(final_prob[i] != NULL);
    }
  

  current = head;
  while(current != NULL)
    {
      /* update_bar(progress_count*100/n_chromosomes); */
      for(int i = 0; i < individual_count; i++)
	{
	  /* update_bar((i+1)*100/individual_count); */
	  for(int j = 0; j < NO_MODELS; j++)
	    {
	      indv_model_prob[i][j] = p_model(current, current->haplotype_1[i], current->haplotype_2[i], j);
	      final_prob[i][j] = final_prob[i][j] + indv_model_prob[i][j];
/* #ifdef DEBUG */
/* 	      printf("%lf\t",indv_model_prob[i][j]); */
/* #endif */
	    }
	  /* printf("\n"); */
	}
      /* printf("\rCalculating the log-likelihood for chromosome %u:100%% Done", current->chrom_id); */
      /* fflush(stdout); */
      current = current->next;
      /* progress_count++; */
      /* printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); */
    }

/*   for(int i = 0; i < individual_count; i++) */
/*     { */
/* /\* #ifdef DEBUG *\/ */
/*       for(int j = 0; j < NO_MODELS; j++) */
/* 	{ */
/* 	  printf("%15lf\t",final_prob[i][j]); */
/* 	} */
/* /\* #endif *\/ */
/*       printf("\n"); */
/*     } */



  char *row_probabilities[individual_count];
  char *new[individual_count];

  for(int i = 0; i < individual_count; i++)
    {
      row_probabilities[i] = malloc(BUFFER_SIZE);
      new[i] = malloc(BUFFER_SIZE);
      assert(row_probabilities[i] != NULL);
      assert(new[i] != NULL);
    }
  


  /* double log_posterior_denominator[individual_count], scale[individual_count]; */
  double *log_posterior_denominator, *scale;
  log_posterior_denominator = malloc((individual_count) * sizeof(double));
  scale = malloc((individual_count) * sizeof(double));
  assert(log_posterior_denominator != NULL);
  assert(scale != NULL);
  
  for(int i = 0; i < individual_count; i++)
    {
      scale[i] = final_prob[i][0];
      for(int j = 1; j < NO_MODELS; j++)
      	{
	  scale[i] = fmax(final_prob[i][j],scale[i]);
	}
      scale[i] = fabs(scale[i])-0.5;
    }



    for(int i = 0; i < individual_count; i++)
    {

      for(int j = 0; j < NO_MODELS; j++)
	{
      log_posterior_denominator[i] = log(exp(final_prob[i][0] + scale[i]) + exp(final_prob[i][1] + scale[i]) + exp(final_prob[i][2] + scale[i]) + exp(final_prob[i][3] + scale[i]) + exp(final_prob[i][4] + scale[i]) + exp(final_prob[i][5] + scale[i])) - scale[i];


	}

    }

    /* printf("The log posterior probabilities ~~~~~~~~~~~~~~~~~~~\n\n"); */
  
  double *log_posterior[individual_count];
  double *posterior_prob[individual_count];
  for(int i = 0; i < individual_count; i++)
    {
      log_posterior[i] = (double*) calloc(NO_MODELS, sizeof(double));
      posterior_prob[i] = (double*) calloc(NO_MODELS, sizeof(double));
      assert(log_posterior[i] != NULL);
      assert(posterior_prob[i] != NULL);
    }

  
  for(int i = 0; i < individual_count; i++)
    {
      for(int j = 0; j < NO_MODELS; j++)
	{
	  log_posterior[i][j] = final_prob[i][j] - log_posterior_denominator[i];
	  posterior_prob[i][j] = exp(log_posterior[i][j]);
	  /* printf("%15lf\t",log_posterior[i][j]); */
	}
      /* printf("\n"); */
    }


  GHashTable *hash_model = g_hash_table_new(g_int_hash, g_int_equal);
  char *model_name_array;
  model_name_array = malloc((NO_MODELS+1) * sizeof(char));
  assert(model_name_array != NULL);
  strcpy(model_name_array,"abcdef");
  int *key;
  key = malloc(NO_MODELS * sizeof(int));
  assert(key != NULL);
  for(int i = 0; i < NO_MODELS; i++)
    {                                   
      key[i] = i;
      g_hash_table_insert(hash_model,&key[i],&model_name_array[i]);
    }                                   
  /* char best_model_post_prob[individual_count+1], second_best_model_post_prob[individual_count+1]; */


  char *best_model_post_prob = NULL;
  best_model_post_prob = malloc(individual_count+1);
  assert(best_model_post_prob != NULL);
 
  char *second_best_model_post_prob = NULL;
  second_best_model_post_prob = malloc(individual_count+1);
  assert(second_best_model_post_prob != NULL);
  
  double *BayesFactor;
  BayesFactor = malloc((individual_count) * sizeof(double));
  assert(BayesFactor != NULL);
  
  for(int i = 0; i < individual_count; i++)
    {
#ifdef DEBUG
      for(int j = 0; j < NO_MODELS; j++)
	{
	  printf("%15lf\t",final_prob[i][j]);
	}
#endif
      largest(posterior_prob[i],NO_MODELS,hash_model,&best_model_post_prob[i],
	      &second_best_model_post_prob[i],&BayesFactor[i]);
      /* printf("\n"); */
    }



  for(int i = 0; i < individual_count; i++)
    {
      for(int j = 0; j < NO_MODELS; j++)

	{
	  /* #ifdef DEBUG */
	  if(j == 0){
	    sprintf(row_probabilities[i],"%15lf\t",final_prob[i][j]);
	  }
	  else{
	    sprintf(new[i],"%15lf\t",final_prob[i][j]);
	    strcat(row_probabilities[i],new[i]);
	  }
	}
      /* strcat(row_probabilities[i],"\n"); */
      /* #endif */
    }


  
  
  /* char append_posterior[individual_count][BUFFER_SIZE]; */
  /* char model_names_BF[individual_count][BUFFER_SIZE]; */

  char *append_posterior[individual_count];
  char *model_names_BF[individual_count];

  for(int i = 0; i < individual_count; i++)
    {
      append_posterior[i] = malloc(BUFFER_SIZE);
      model_names_BF[i] = malloc(BUFFER_SIZE);
      assert(append_posterior[i] != NULL);
      assert(model_names_BF[i] != NULL);
    }

  double *entropy;
  entropy = calloc((individual_count), sizeof(double));
  assert(entropy != NULL);
  for(int i = 0; i < individual_count; i++)
    {
      for(int j = 0; j < NO_MODELS; j++)
	{
	  if(posterior_prob[i][j] != 0)
	    {
	      entropy[i] = entropy[i] - (posterior_prob[i][j] * log2(posterior_prob[i][j]));
	    }
	}
      /* printf("%15lf\n",entropy[i]); */
    }

  
  for(int i = 0; i < individual_count; i++)
    {
      for(int j = 0; j < NO_MODELS; j++)
  	{
	  sprintf(append_posterior[i],"%15lf\t",posterior_prob[i][j]);
	  strcat(row_probabilities[i],append_posterior[i]);
  	}
      sprintf(model_names_BF[i],"\t%8c\t%8c\t%15lf\t%15lf",best_model_post_prob[i],
	      second_best_model_post_prob[i], BayesFactor[i], entropy[i]);
      strcat(row_probabilities[i],model_names_BF[i]);
      strcat(row_probabilities[i],"\n");
    }


  char header[BUFFER_SIZE] = "log_prob_A";
  strcat(header,"\tlog_prob_B");
  strcat(header,"\tlog_prob_C");
  strcat(header,"\tlog_prob_D");
  strcat(header,"\tlog_prob_E");
  strcat(header,"\tlog_prob_F");

  strcat(header,"\tpost_prob_A");
  strcat(header,"\tpost_prob_B");
  strcat(header,"\tpost_prob_C");
  strcat(header,"\tpost_prob_D");
  strcat(header,"\tpost_prob_E");
  strcat(header,"\tpost_prob_F");
  strcat(header,"\tbest_model");
  strcat(header,"\tnext_best");
  strcat(header,"\tBayes_Factor");
  strcat(header,"\tEntropy");
  strcat(header,"\n");
  /* printf("%s\n", header); */
  
  gsize bytes_written;
  /* outfile = g_io_channel_new_file(argv[9],"w",&error); */
  outfile = g_io_channel_new_file(output_filename,"w",&error);
  g_io_channel_write_chars(outfile,header,strlen(header),&bytes_written,&error);
  for(int i = 0; i < individual_count; i++)
    {

      g_io_channel_write_chars(outfile,row_probabilities[i],strlen(row_probabilities[i]),&bytes_written,&error);
    }





  
  
  for(int i = 0; i < individual_count; i++)
    {
      free(indv_model_prob[i]);
      free(final_prob[i]);
      free(log_posterior[i]);
      free(posterior_prob[i]);

    }


  free(log_posterior_denominator);
  free(scale);


  free(model_name_array);
  free(key);                                                                                    
  g_hash_table_destroy(hash_model);


  for(int i = 0; i < individual_count; i++)
    {
      free(row_probabilities[i]);
      free(new[i]);
      free(append_posterior[i]);
      free(model_names_BF[i]);
    }

  free(entropy);
  free(BayesFactor);
  free(best_model_post_prob);
  free(second_best_model_post_prob);
  
  
  delete_chrom_data(head);
  free(chrom_index);
  free(chrom_no_loci);
  free(chrom_length);
  free(chrom_recom_rate);
  free(loci);


  g_io_channel_unref(infile_chrom);
  g_io_channel_unref(infile_A);
  g_io_channel_unref(infile_B);
  g_io_channel_unref(infile_indv);

  g_io_channel_unref(outfile);
  
  return 0;
  
}
