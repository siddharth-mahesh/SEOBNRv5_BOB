#include"sBOB_funcs.h"

#define MAXCHAR 1000

int main(){
	double t,hp,hc,amp,phase;
	double Mf, af;

	af = 0.6864;
	Mf = 0.9516;

	double *h , *ph;
	const double t_init = -100.0;
	const double t_end = 100.0;
	const double delta_t = 0.1;
	FILE *outfile;
	char str[MAXCHAR];
	
	
	t = t_init;
	
	outfile = fopen("sBOB_waveform.txt","w+");
	h = &amp;
	ph = &phase;
	do{
		//printf("t");
		get_sBOB_strainamplitude(t,Mf,af,h);
		get_sBOB_phase(t,Mf,af,ph);
		amp = *h;
		phase = *ph;
		hp = amp*cos(phase);
		hc = amp*sin(phase);
		fprintf(outfile,"%f \t %f \t %f \t %f \t %f \n",t,hp,hc,amp,phase);
		t = t + delta_t;
	}while (t < t_end);
	
	fclose(outfile);
	return 0;	
}