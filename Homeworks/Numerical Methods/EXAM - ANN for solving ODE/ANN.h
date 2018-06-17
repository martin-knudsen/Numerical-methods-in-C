#ifndef HAVE_ANN
#define HAVE_ANN
typedef struct {int n; double (*f)(double); double (*fm)(double); gsl_vector* data;} ann;
ann* ann_alloc(int number_of_hidden_neurons, double(*activation_function)(double), double(*activation_function_der)(double));
void ann_free(ann* network);
void ann_feed_forward(ann* network, double* F, double* dF);
void ann_train(ann* network, gsl_vector* xlist,double (*function)(double) ,double x0, double y0);

#endif
