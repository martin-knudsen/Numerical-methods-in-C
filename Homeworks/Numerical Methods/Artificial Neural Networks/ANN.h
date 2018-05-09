#ifndef HAVE_ANN
#define HAVE_ANN
typedef struct {int n; double (*f)(double); gsl_vector* data;} ann;
ann* ann_alloc(int number_of_hidden_neurons, double(*activation_function)(double));
void ann_free(ann* network);
double ann_feed_forward(ann* network, double x);
void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist);

typedef struct {int n; double (*f)(double); gsl_vector* data;} ann2D;
ann2D* ann2D_alloc(int number_of_hidden_neurons, double(*activation_function)(double));
void ann2D_free(ann2D* network);
double ann2D_feed_forward(ann2D* network, double x1, double x2);
void ann2D_train(ann2D* network, gsl_matrix* xlist, gsl_vector* ylist);

#endif
