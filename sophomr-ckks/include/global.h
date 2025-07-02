
int payload_size = 0;// 612 * 8;

int num_transaction = 65536;
int num_pertinent = 50;

int ptxt_modulus = 786433;
int degree = 32768; // 65536도 너무 작다고 이거 2배는 되야한다는데

int NumLargeDigits = 2; 
int ScalingModSize = 50; // 60크다고 오류뜸 // 근데 이것도 올려도 될듯
int MultiplicativeDepth = 38;

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
//param PSparam(1024, ptxt_modulus, 2, 80, 0.5); 

std::vector<std::vector<uint32_t>> coeff_rangeCheck = 
{{0,216615,237416,741379,99137,105134,235513,643136},
{498247,88199,302058,515923,339325,161996,88430,643211},
{342302,233360,731664,655413,219378,239548,347203,445370},
{361791,80486,126394,391207,551639,528197,478787,88772},
{710790,716155,456692,301062,159136,300549,8286,765893},
{1,0,0,0,0,0,0,0}};