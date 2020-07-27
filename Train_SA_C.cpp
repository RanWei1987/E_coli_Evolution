#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;
using namespace RcppEigen;

//[[Rcpp::export]]
NumericMatrix MatrixMul(NumericMatrix x, NumericMatrix y);
NumericMatrix MatrixMul_schur(NumericMatrix x, NumericMatrix y);  
List Train_SA_C(NumericMatrix G, NumericMatrix E, NumericMatrix X_tr,NumericMatrix Y_tr, NumericVector GR_tr, double obj_best, NumericVector T_range, int Iter);

//[[Rcpp::export]]
List Train_SA_C(NumericMatrix G, NumericMatrix E, NumericMatrix X_tr,NumericMatrix Y_tr, NumericVector GR_tr, double obj_best, NumericVector T_range, int Iter){
	cout<<"Firstone: "<<R::runif(0,1)<<endl;

	Environment global = Environment::global_env();
    	int num_cycle = global["num_cycle"];
  	int num_round = global["num_round"];
	cout<<"Cycle "<<num_cycle<<"Round "<<num_round<<endl;
 	
	ofstream outfile;
	outfile.open("run_report.txt", ios::out | ios::app);
	
	outfile<<"Cycle "<<num_cycle<<"Round "<<num_round<<endl;
   	outfile.close();

	int num_exp = GR_tr.length();
	double obj = 1e30, obj_prev = 1, obj_curr = 1;
	int traits = G.nrow(), genes = G.ncol(), envi = E.ncol();
        int numCirc = 1000;
        NumericVector circBuf;
	circBuf = rep(1, numCirc); 	
	double sumCirc = sum(circBuf);
	int iCirc = 0;
        double pBad = 1;

        double SA_Tinit = T_range[0];
        double SA_Tfinal = T_range[1];

//	int ran_tracer = 1;

	double SA_numIters = (double)Iter;
	int SA_iter = 0;
	double SA_s, SA_pAcc, SA_r;
	double SA_T = SA_Tinit;
	double SA_lambda = -log(SA_Tfinal / SA_Tinit);
	List L = List::create(G, E);
	int success = 0, failure = 0;
	NumericVector GR_guess_tr (GR_tr.length(), 0);
	NumericMatrix TG_tr = MatrixMul(G, X_tr);
	NumericMatrix TE_tr = MatrixMul(E, Y_tr);

	float part = (float) genes / (genes + envi);
	
	NumericMatrix G_best = G;
	NumericMatrix E_best = E;
	NumericVector best = NumericVector::create(Named("iter",0), Named("obj")=0, Named("rho")=0, Named("GRdif")=0);
	char select;
	int whichRow, whichCol;
	double delta, sum_x, sum_y, sum_xy, sum_xx, sum_yy, pearson, GR_mean_diff, obj_new;
	int status_freq = 10000;

	while (true){
		float current = (float) R::runif(0, 1);
		whichRow = int (floor(R::runif(0, traits)));
  		if (current <= part){
			select = 'G';
			whichCol = int (floor(R::runif(0, genes))); 
		} else {
			select = 'E';			
			whichCol = int (floor(R::runif(0, envi)));	
		}
		delta = R::runif(-pBad/10, pBad/10);
	//	ran_tracer = ran_tracer+4;
	//	cout<<current<<"  "<<whichRow<<"  "<<whichCol<<"  "<<delta<<endl;
		//if (pBad > 1e-11) {
		//	delta = R::runif(-1, 1);
		//}
	 //	cout<<current<<'\t'<<part<<endl;
		if (select =='G'){
			G( whichRow, whichCol ) += delta;
     			TG_tr( whichRow, _ ) = TG_tr( whichRow, _ ) + delta * X_tr( whichCol, _ );
		} else {
			E( whichRow, whichCol ) += delta;
     			TE_tr( whichRow, _ ) = TE_tr( whichRow, _ ) + delta * Y_tr( whichCol, _ );
    		}
		NumericMatrix product = MatrixMul_schur(TG_tr, TE_tr);

 		for (int i=0; i<num_exp; i++) {
	         	GR_guess_tr[i] = sum(product( _, i )) / 1500;
		}

		sum_x = sum(GR_tr);
		sum_y = sum(GR_guess_tr);
		sum_xy = sum(GR_tr * GR_guess_tr);
		sum_xx = sum(pow(GR_tr, 2));
		sum_yy = sum(pow(GR_guess_tr, 2));
 		pearson = (num_exp * sum_xy - sum_x * sum_y) / sqrt( (num_exp * sum_xx - sum_x * sum_x) * (num_exp * sum_yy - sum_y * sum_y) );
		
		GR_mean_diff = mean(abs(GR_guess_tr - GR_tr));
		
		obj_new = abs(GR_mean_diff * (1-pearson));
		
		SA_s = SA_iter / SA_numIters;
		SA_iter++;
		SA_T = SA_Tinit * exp (-SA_lambda * SA_s);
		SA_pAcc = exp((obj-obj_new) / SA_T);
		SA_r = R::runif(0, 1);
//		ran_tracer++;
		// cout<<SA_T<<endl;
		//if (SA_iter % 10000 == 0){
		//	cout<<ran_tracer<<'\t'<<endl;
	        //}	
		// cout<<SA_s<<'\t'<<SA_iter<<'\t'<<SA_numIters<<endl;
		//cout<<obj<<'\t'<<obj_new<<'\t'<<((obj-obj_new) / SA_T)<<'\t'<<SA_pAcc<<endl;
	        // cout<<SA_r<<"   "<<SA_pAcc<<obj_new<<SA_T<<endl;	
		if (SA_r < SA_pAcc) {
	  		if (obj > obj_new){
				G_best = G;
        			E_best = E;
       				best["iter"] = SA_iter;
        			best["obj"] = obj_new;
        			best["rho"] = pearson;
        			best["GRdif"] = GR_mean_diff;
				success++;
		      	}
			obj = obj_new;
		} else {
			if (select =='G'){
				G( whichRow, whichCol ) -= delta;
     				TG_tr( whichRow, _ ) = TG_tr( whichRow, _ ) - delta * X_tr( whichCol, _ );
			} else {
				E( whichRow, whichCol ) -= delta;
     				TE_tr( whichRow, _ ) = TE_tr( whichRow, _ ) - delta * Y_tr( whichCol, _ );
    			}	
			iCirc++;
			sumCirc -= circBuf[iCirc % numCirc];
			circBuf[iCirc % numCirc] = SA_pAcc;
			sumCirc += SA_pAcc;
			failure++;
		        pBad = sumCirc / numCirc;
		}
	//	if (abs(obj_curr - obj_prev) > 0.5){
	//		status_freq = 100;
	//	} else {
	//		status_freq = 10000;
	//	}
//		if (ran_tracer == 96666){
//			cout<<ran_tracer<<'\t'<<R::runif(0,1)<<endl;
//		}
		if (SA_iter % status_freq == 0) {
			cout<<"Iter "<<SA_iter<<" pBad "<<pBad<<" T "<<SA_T<<" Obj "<<obj<<" Rho "<<pearson<<" GRdif "<<GR_mean_diff<<" success "<<success<<" failure "<<failure<<" obj_new "<<obj_new<<" SA_pAcc "<<SA_pAcc<<" Delta "<<delta<<" SA_r "<<SA_r<<endl;
//			cout<<"tracer: "<<ran_tracer<<'\t'<<"random value: "<<R::runif(0,1)<<endl;
 			ofstream outfile;
			outfile.open("run_report.txt", ios::out | ios::app);
			outfile<<"Iter "<<SA_iter<<" pBad "<<pBad<<" T "<<SA_T<<" Obj "<<obj<<" Rho "<<pearson<<" GRdif "<<GR_mean_diff<<" success "<<success<<" failure "<<failure<<" obj_new "<<obj_new<<" SA_pAcc "<<SA_pAcc<<" Delta "<<delta<<" SA_r "<<SA_r<<endl;
   			outfile.close();
			
			if (SA_iter >= SA_numIters) {
				cout<<"Cannot move further"<<endl;
	 			cout<<best<<endl;
				NumericMatrix pred_current = MatrixMul_schur(MatrixMul(G_best, X_tr), MatrixMul(E_best, Y_tr)) / 1500;
 				NumericVector pred_rate(pred_current.ncol());
				for (int i = 0; i < pred_current.ncol(); i++) {
					pred_rate[i] = sum(pred_current( _, i ));
				}
				ofstream outfile;
				outfile.open("run_report.txt", ios::out | ios::app);
				outfile<<"Cannot move further"<<endl;
				outfile<<"iter "<<best[0]<<" obj"<<best[1]<<" rho "<<best[2]<<" GRdif "<<best[3]<<endl;
				outfile<<"Lab Growth Rates:"<<'\n'<<GR_tr<<endl;
				outfile<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<pred_rate<<endl;
   				outfile.close();
				
				cout<<best<<endl;
				cout<<"Train: "<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_tr<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<GR_guess_tr<<endl;	
				
				num_cycle++;
				global["num_cycle"] = num_cycle;
				List runagain;
				if (num_cycle>3) {
					num_cycle = 1;
					global["num_cycle"] = 1;
					num_round++;

					if (num_round > 3) {
						cout<<"Cannot get a good solution. Raise the goal of obj by 100%"<<endl;
						obj_best = obj_best * 2;
						num_round = 1;
						global["num_round"] = 1;
						cout<<"Now the goal of obj is: "<<obj_best<<'\n'<<"Round "<<num_round<<"Cycle "<<num_cycle<<endl;

						ofstream outfile;
						outfile.open("run_report.txt", ios::out | ios::app);
						outfile<<"Cannot get a good solution. Raise the goal of obj by 100%"<<'\n';
						outfile<<"Now the goal of obj is: "<<obj_best<<'\n'<<"Round "<<num_round<<"Cycle "<<num_cycle<<endl;
   						outfile.close();
					}

						cout<<"G, E redefined. Anneal again with new G and E as input"<<endl;
						
						ofstream outfile;
		 				outfile.open("run_report.txt", ios::out | ios::app);
						outfile<<"G, E redefined. Anneal again with new G and E as input"<<endl;
   						outfile.close();
					
						NumericVector G_initial = rnorm(traits * genes);
						G_initial.attr("dim") = Dimension(traits, genes);
						global["G_initial"] = G_initial;
						
						NumericVector E_initial = rnorm(traits * envi);
						E_initial.attr("dim") = Dimension(traits, envi);
						global["E_initial"] = E_initial;

						G = as<NumericMatrix>(G_initial);	
						E = as<NumericMatrix>(E_initial);
			
						ofstream outfile_2;
		 				outfile_2.open("run_result.txt", ios::out | ios::app);
						outfile_2<<"Initial G: "<<G<<'\n'<<"Initial E: "<<E<<endl;
   						outfile_2.close();

 						runagain = Train_SA_C(G, E, X_tr, Y_tr, GR_tr, obj_best, T_range, Iter);
				} else{
					cout<<"Anneal again with the best G and E so far as input"<<endl;							
					ofstream outfile;
		 			outfile.open("run_report.txt", ios::out | ios::app);
					outfile<<"Anneal again with the best G and E so far as input"<<endl;
   					outfile.close();

					runagain = Train_SA_C(G, E, X_tr, Y_tr, GR_tr, obj_best, T_range, Iter);
				}
				return(runagain);	
			}
		}
//		if (iter == 5870000) {
//			status_freq = 100;
//		}	
		if (obj < obj_best) {
			cout<<"Good enough"<<'\n';
			cout<<"iter "<<SA_iter<<"obj "<<obj<<"rho "<<pearson<<"GRdif "<<GR_mean_diff<<'\n';
			cout<<"Train: "<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_tr<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<GR_guess_tr<<endl;	
		
			ofstream outfile;
			outfile.open("run_report.txt", ios::out | ios::app);
			outfile<<"Good enough"<<'\n';
			outfile<<"iter "<<SA_iter<<"obj "<<obj<<"rho "<<pearson<<"GRdif "<<GR_mean_diff<<'\n';

			outfile<<"Lab Growth Rates:"<<'\n'<<GR_tr<<'\n';
			outfile<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<GR_mean_diff<<endl;
	   		outfile.close();
			return (List::create(obj_new, G, E, L, GR_guess_tr, SA_iter, GR_mean_diff));		
		}				
	}
}

// [[Rcpp::depends(RcppEigen)]]
NumericMatrix MatrixMul(NumericMatrix x, NumericMatrix y) {
  Eigen::Map<Eigen::MatrixXd> X = as<Eigen::Map<Eigen::MatrixXd> >(x);
  Eigen::Map<Eigen::MatrixXd> Y = as<Eigen::Map<Eigen::MatrixXd> >(y);
  Eigen::MatrixXd Z = X * Y;
  return (Rcpp::NumericMatrix(wrap(Z)));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix MatrixMul_schur(NumericMatrix x, NumericMatrix y) { 
  arma::mat X = as<arma::mat>(x);
  arma::mat Y = as<arma::mat>(y); 
  arma::mat Z = X % Y;
  return(Rcpp::NumericMatrix(wrap(Z)));
} 
