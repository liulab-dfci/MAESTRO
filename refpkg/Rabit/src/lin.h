/*
 * lin.h
 *
 *  Created on: Nov 2, 2014
 *      Author: peng
 */
#ifndef LIN_H_
#define LIN_H_


// X is column major order
void QR_condense(double X[], double Y[], const int n, const int p, double *Y_last);


#endif /* LIN_H_ */
