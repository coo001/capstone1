
int payload_size = 0;// 612 * 8;

int num_transaction = 2048;
int num_pertinent = 50;

int ptxt_modulus = 786433;
int degree = 1024; // 65536�� �ʹ� �۴ٰ� �̰� 2��� �Ǿ��Ѵٴµ�

int NumLargeDigits = 3; 
int ScalingModSize = 50;// 50; // 60ũ�ٰ� ������ // �ٵ� �̰͵� �÷��� �ɵ�
int MultiplicativeDepth = 60;

int NumLargeDigits_comp = 1;
int MultiplicativeDepth_comp = 12;
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


int range = 39;