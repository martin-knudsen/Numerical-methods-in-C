#ifndef HAVE_ANN
#define HAVE_ANN
typedef struct {int n; double (*f)(double); gsl_vector* data;} ann;
ann* ann_alloc(int number_of_hidden_neurons, double(*activation_function)(double));
void ann_free(ann* network);
double ann_feed_forward(ann* network, double x);
void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist);

#endif
