
int payload_size = 0;//612 * 8;

int num_transaction = 65536;
int num_pertinent = 20;

int ptxt_modulus = 786433;
int degree = 65536;

int NumLargeDigits = 2;
int ScalingModSize = 60;
int MultiplicativeDepth = 40; 
// r이 이거 때매 안 됐네

int NumLargeDigits_comp = 1;
int MultiplicativeDepth_comp = 1;
int dim_trace = 4;

int payload_len;
int numctxt;
int degree_half;
int b_tilde1, g_tilde1;
int numrow, numrow_po2;
int b_tilde2, g_tilde2;
int degree_trace;
int degree_trace_half;
int trace_swap;
int trace_shift;

param PSparam(256, ptxt_modulus, 2, 80, 0.5); 


int range_r = 300;
std::vector<std::vector<uint32_t>> coeff_rangeCheck; 
/*=
{ 
{0,232953,198154,758413,195212,181399,772174,658710},
{104565,151354,451694,579908,234693,668615,231936,391627},
{317379,521977,569697,386053,399731,721882,722879,263501},
{394458,559203,749800,195862,146633,65481,454739,1527},
{6890,695671,605212,54488,487754,577732,412310,628533},
{764293,1,0,0,0,0,0,0}
};*/
