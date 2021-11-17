package bartMachine;

import java.io.Serializable;
import org.apache.commons.math.special.Gamma;

import gnu.trove.list.array.TIntArrayList;

/**
 * This portion of the code implements the informed prior information on covariates feature.
 * 
 * @author 	Adam Kapelner and Justin Bleich
 * @see 	Section 4.10 of Kapelner, A and Bleich, J. bartMachine: A Powerful Tool for Machine Learning in R. ArXiv e-prints, 2013
 */
public class bartMachine_i_prior_cov_spec extends bartMachine_h_eval implements Serializable{
	
	/** Do we use this feature in this BART model? */
	protected boolean use_prior_cov_spec;
	/** This is a probability vector which is the prior on which covariates to split instead of the uniform discrete distribution by default */
	protected double[] cov_split_prior;	



	/**
	 * Pick one predictor from a set of valid predictors that can be part of a split rule at a node
	 * while accounting for the covariate prior.
	 * 
	 * @param node	The node of interest
	 * @return		The index of the column to split on
	 */
	private int pickRandomPredictorThatCanBeAssignedF1(bartMachineTreeNode node){
		TIntArrayList predictors = node.predictorsThatCouldBeUsedToSplitAtNode();
		//get probs of split prior based on predictors that can be used and weight it accordingly
		double[] weighted_cov_split_prior_subset = getWeightedCovSplitPriorSubset(predictors);
		//choose predictor based on random prior value	
		return StatToolbox.multinomial_sample(predictors, weighted_cov_split_prior_subset);
	}
	
	/**
	 * The prior-adjusted number of covariates available to be split at this node
	 *  
	 * @param node		The node of interest
	 * @return			The prior-adjusted number of covariates that can be split
	 */
	private double pAdjF1(bartMachineTreeNode node) {
		if (node.padj == null){
			node.padj = node.predictorsThatCouldBeUsedToSplitAtNode().size();
		}
		if (node.padj == 0){
			return 0;
		}
		if (node.isLeaf){
			return node.padj;
		}			
		//pull out weighted cov split prior subset vector
		TIntArrayList predictors = node.predictorsThatCouldBeUsedToSplitAtNode();
		//get probs of split prior based on predictors that can be used and weight it accordingly
		double[] weighted_cov_split_prior_subset = getWeightedCovSplitPriorSubset(predictors);	
		
		//find index inside predictor vector
		int index = bartMachineTreeNode.BAD_FLAG_int;
		for (int i = 0; i < predictors.size(); i++){
			if (predictors.get(i) == node.splitAttributeM){
				index = i;
				break;
			}
		}
		
		//return inverse probability		
		return 1 / weighted_cov_split_prior_subset[index];
	}
	
	/**
	 * Given a set of valid predictors return the probability vector that corresponds to the
	 * elements of <code>cov_split_prior</code> re-normalized because some entries may be deleted
	 * 
	 * @param predictors	The indices of the valid covariates
	 * @return				The updated and renormalized prior probability vector on the covariates to split
	 */
	private double[] getWeightedCovSplitPriorSubset(TIntArrayList predictors) {
		double[] weighted_cov_split_prior_subset = new double[predictors.size()];
		for (int i = 0; i < predictors.size(); i++){
			weighted_cov_split_prior_subset[i] = cov_split_prior[predictors.get(i)];
		}
		Tools.normalize_array(weighted_cov_split_prior_subset);
		return weighted_cov_split_prior_subset;
	}	

	public void setCovSplitPrior(double[] cov_split_prior) {
		this.cov_split_prior = cov_split_prior;
		//if we're setting the vector, we're using this feature
		use_prior_cov_spec = true;
	}

	
	/////////////nothing but scaffold code below, do not alter!
	
	public int pickRandomPredictorThatCanBeAssigned(bartMachineTreeNode node){
		if (use_prior_cov_spec){
			return pickRandomPredictorThatCanBeAssignedF1(node);
		}
		return super.pickRandomPredictorThatCanBeAssigned(node);
	}	
	
	public double pAdj(bartMachineTreeNode node){
		if (use_prior_cov_spec){
			return pAdjF1(node);
		}
		return super.pAdj(node);
	}

    // protected void updateCovSplit(int sample_num, bartMachineTreeNode[] bart_trees) {
        
    // }


    // TONY
    public void initGibbsSamplesCovSplit(double[][] gibbs_samples_cov_split_prior) {
        this.gibbs_samples_cov_split_prior = gibbs_samples_cov_split_prior; 
    }

    protected void updateCovSplit(int sample_num, bartMachineTreeNode[] bart_trees) {
        // System.out.println("do_ard = " + do_ard);
        double b = 5.0 * (double)p / alpha_0;
        int num_anneal = num_gibbs_burn_in;
        int num_anneal_top = (int)Math.floor((double)num_gibbs_burn_in / 2.0);
        double a = 1.0 / ((double)num_anneal-1) * (1.0 - b);

        if(do_ard == 1 && sample_num < num_anneal && sample_num > num_anneal_top) {

            if(do_prior == 0) {
                double t_star = a * sample_num + b; 
                double alpha_star = t_star * alpha_0;
                updateAlpha(bart_trees, alpha_star);
                cov_split_prior = StatToolbox.sample_from_dirichlet(alpha_up);
            }
            
            else {
                updateAlpha(bart_trees, alpha_0);
                log_cov_split_prior = sample_log_dirichlet(alpha_up);
                cov_split_prior = Exp(log_cov_split_prior);
                updateAlpha0();
            }

        }

        if(do_ard == 1 && sample_num >= num_anneal) {

            if(do_prior == 0) {
                updateAlpha(bart_trees, alpha_0);
                cov_split_prior = StatToolbox.sample_from_dirichlet(alpha_up);
            }

            else {
                updateAlpha(bart_trees, alpha_0);
                log_cov_split_prior = sample_log_dirichlet(alpha_up);
                cov_split_prior = Exp(log_cov_split_prior);
                updateAlpha0();
                // if(sample_num % 100 == 0) System.gc();
            }
        }

        gibbs_samples_cov_split_prior[sample_num-1] = cov_split_prior.clone();
        // if(sample_num % 100 == 0) System.gc();
        // if(sample_num % 100 == 0)
            // System.out.println("alpha_0 = " + alpha_0);
    }

    protected double[] Exp(double[] x) {
        double[] out = new double[x.length];
        for(int i = 0; i < x.length; i++)
            out[i] = Math.exp(x[i]);
        return(out);
    }

    protected double[] sample_log_dirichlet(double[] alpha) {
        double[] out = new double[alpha.length];

        for(int i = 0; i < alpha.length; i++) {
            out[i] = StatToolbox.sample_from_gamma(alpha[i], 1.0);
            while(out[i] == 0) {
                out[i] = StatToolbox.sample_from_gamma(alpha[i], 1.0);
            }
            out[i] = Math.log(out[i]);
        }

        double log_sum_exp = LogSumExp(out);
        for(int i = 0; i < alpha.length; i++) {
            out[i] -= log_sum_exp;
        }

        return(out);
    }

    // TONY : Define updateAlpha() 
    // TONY : To implement this, I need to get the variable counts.
    protected void updateAlpha(bartMachineTreeNode[] bart_trees, double alpha_0_) {
        int[] counts = new int[p]; 
        for(bartMachineTreeNode tree : bart_trees) {
            counts = Tools.add_arrays(counts, tree.attributeSplitCounts2());
        }
        for(int i = 0; i < p; i++) {
            alpha_up[i] = alpha_0_ / (double)p + (double)counts[i]; // TONY TODO : Replace 1.0 with something else
        }
        counts = null;
    } 


    // TONY : Below here is prior on alpha

    protected double BetaPrior(double alpha) {
        double u = alpha / (alpha + (double)p);
        // double u = alpha / (alpha + 1.0);
        double a = 0.5;
        // double a = 1.0;
        // double a = 1.0 / (double)p;
        double b = 1.0;

        return (a - 1.0) * Math.log(u) + (b - 1.0) * Math.log(1-u);
        // return (a - 1.0) * Math.log(alpha) - (a + b) * Math.log(alpha);
    }

    protected double loglik_cov_prior = 0.0;
    protected void CalcSumLogCovPrior() {
        loglik_cov_prior = 0.0;
        for(int i = 0; i < p; i++) {
            loglik_cov_prior += log_cov_split_prior[i];
        }
        loglik_cov_prior /= p;
    }

    protected double AlphaLoglik(double alpha) {
        // System.out.println("alpha = " + alpha);
        // System.out.println("logGamma(alpha) = " + Gamma.logGamma(alpha));
        // System.out.println("logGamma(alpha/p) = " + Gamma.logGamma(alpha / (double)p));
        // System.out.println("loglik_cov_prior = " + loglik_cov_prior);
        return alpha * loglik_cov_prior + Gamma.logGamma(alpha) - Gamma.logGamma(alpha / (double)p) * (double)p;
    }

    protected double Max(double[] x) {
        double retval = x[0];
        for(int i = 1; i < x.length; i++) {
            if(x[i] > retval)
                retval = x[i];
        }
        return retval;
    }


    protected double LogSumExp(double[] x) {
        double M = Max(x);
        double sumval = 0.0;

        for(int i = 0; i < x.length; i++)
            sumval += Math.exp(x[i] - M);

        return M + Math.log(sumval);
    }

    protected int SampleDiscrete(double[] w) {
        double u = StatToolbox.rand();
        // System.out.println("U = " + u);
        double cs = 0.0;
        int out = w.length;

        for(int i = 0; i < w.length; i++) {
            cs += w[i];
            // System.out.println("cs = " + cs);
            if(u < cs) {
                out = i; 
                break;
            }
        }
        return out;
    }

    int NUM_CAND = 1000;
    protected void updateAlpha0() {
        // First get the candidate values, discrete grid, etc
        // System.out.println("Fill in candidates");
        double[] candidates = new double[NUM_CAND];
        double[] log_weight = new double[NUM_CAND];
        for(int i = 0; i < NUM_CAND; i++) {
            candidates[i] = (double)(i + 1) / (double)(NUM_CAND+1);
            candidates[i] = (double)p * candidates[i] / (1.0 - candidates[i]);
            // candidates[i] = candidates[i] / (1.0 - candidates[i]);
        }

        // System.out.println("Calc loglik");
        // Update the loglik
        CalcSumLogCovPrior();

        // Fill in log_weights
        // System.out.println("Fill in weights");
        for(int i = 0; i < NUM_CAND; i++) {
            log_weight[i] = BetaPrior(candidates[i]) + AlphaLoglik(candidates[i]);
            // System.out.println("BetaPrior[" + i + "]: " + BetaPrior(candidates[i]));
            // System.out.println("AlphaLoglik[" + i + "]: " + AlphaLoglik(candidates[i]));
        }

        // Get Noramlization and normalize
        // System.out.println("Get normal");
        double log_sum_exp = LogSumExp(log_weight);
        for(int i = 0; i < NUM_CAND; i++) {
            log_weight[i] = Math.exp(log_weight[i] - log_sum_exp);
            // System.out.println(log_weight[i]);
        }
        
        // Sample the value of alpha
        // System.out.println("Sample");

        int update_alpha = SampleDiscrete(log_weight);
        alpha_0 = candidates[update_alpha];

        // System.out.println("set value");
        // System.out.println("updated value = " + update_alpha);
    }
}
