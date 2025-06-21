
int payload_size = 0;// 612 * 8;

int num_transaction = 65536;
int num_pertinent = 10;
int num_person = 4;

int ptxt_modulus = 786433;
int degree = 65536;

int NumLargeDigits = 2;
int ScalingModSize = 60;

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

std::vector<int> num_set;
int id_num;
int range_r;
int MultiplicativeDepth = 40; // 동형 연산 깊이
std::vector<std::vector<uint32_t>> coeff_rangeCheck;