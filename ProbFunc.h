//
// Created by Nan on 7/24/2017.
//

#ifndef OVERLAPFIND_PROBFUNC_H
#define OVERLAPFIND_PROBFUNC_H

#endif //OVERLAPFIND_PROBFUNC_H

long int statistical_bound_of_waiting_time1(double p, long int k, double alpha);
long int statistical_bound_of_waiting_time2(double p, long int k, double alpha);
double *randomwalk_probability_of_pos1(double pI, long int L);
double *randomwalk_probability_of_pos2(double pI, long int L);
double *randomwalk_probability_of_pos3(double pI, long int L);
long int statistical_bound_of_randomwalk1(double pI, long int L, double alpha);
long int statistical_bound_of_randomwalk2(double pI, long int L, double alpha);
double Evalue(double K, double Lambda, long int query_size, long int text_size,
              long int score);
long int MinScore(double K, double Lambda, long int query_size, long int text_size,
                  double MaxEvalue);
double BitScore(double K, double Lambda, long int score);
double P_mutation_bias(long int freq[]);

double ** computeLettersFrequency(long int nb_letters[2][4]);
double ** computeBackgroundFrequency(double ** letters_frequency /* [2][4] */);
double ** computeBackgroundTripletFrequency(long int nb_words[2][64]);

double computeLambda(double ** freq_letters /* [2][4] */);
double computeK(double ** freq_letters /* [2][4] */, double lambda);

double entropyTriplet(long int * count /* [4^3]*/);
double mutualInformationTriplet(long int ** count /* [4^3] x [4^3] */);

