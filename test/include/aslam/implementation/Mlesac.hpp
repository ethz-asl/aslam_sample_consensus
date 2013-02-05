
namespace aslam {

    template<typename P>
    Mlesac<P>::Mlesac(int maxIterations, double threshold, double probability) :
        SampleConsensus<P>(maxIterations, threshold, probability),
        iterations_EM_(3)
    {

    }

    template<typename P>
    Mlesac<P>::~Mlesac(){}


    template<typename PROBLEM_T>
    bool Mlesac<PROBLEM_T>::computeModel(int debug_verbosity_level)
    {
        // Warn and exit if no threshold was set
        if (threshold_ == std::numeric_limits<double>::max())
        {
            SM_ERROR ("[sm::MaximumLikelihoodSampleConsensus::computeModel] No threshold set!\n");
            return (false);
        }

        iterations_ = 0;
        double d_best_penalty = std::numeric_limits<double>::max();
        double k = 1.0;

        std::vector<int> best_model;
        std::vector<int> selection;
        model_t model_coefficients;
        std::vector<double> distances;

        // Compute sigma - remember to set threshold_ correctly !
        //sigma_ = computeMedianAbsoluteDeviation (sac_model_->getInputCloud (), sac_model_->getIndices (), threshold_);
        /// \todo figure out what is sigma
        double sigma_ = 0.1;
        if (debug_verbosity_level > 1)
            SM_DEBUG ("[sm::MaximumLikelihoodSampleConsensus::computeModel] Estimated sigma value: %f.\n", sigma_);

        // Compute the bounding box diagonal: V = sqrt (sum (max(pointCloud) - min(pointCloud)^2))
        //Eigen::Vector4f min_pt, max_pt;
        //getMinMax (sac_model_->getInputCloud (), sac_model_->getIndices (), min_pt, max_pt);
        //max_pt -= min_pt;
        //double v = sqrt (max_pt.dot (max_pt));
        // \todo Figure out how to set V
        double v = 1.0;

        int n_inliers_count = 0;
        size_t indices_size;
        unsigned skipped_count = 0;
        // supress infinite loops by just allowing 10 x maximum allowed iterations for invalid model parameters!
        const unsigned max_skip = max_iterations_ * 10;
  
        // Iterate
        while (iterations_ < k && skipped_count < max_skip)
        {
            // Get X samples which satisfy the model criteria
            sac_model_->getSamples (iterations_, selection);

            if (selection.empty ()) break;

            // Search for inliers in the point cloud for the current plane model M
            if (!sac_model_->computeModelCoefficients (selection, model_coefficients))
            {
                //iterations_++;
                ++ skipped_count;
                continue;
            }

            // Iterate through the 3d points and calculate the distances from them to the model
            sac_model_->getDistancesToModel (model_coefficients, distances);

            // Use Expectiation-Maximization to find out the right value for d_cur_penalty
            // ---[ Initial estimate for the gamma mixing parameter = 1/2
            double gamma = 0.5;
            double p_outlier_prob = 0;

            indices_size = sac_model_->getIndices ()->size ();
            std::vector<double> p_inlier_prob (indices_size);
            for (int j = 0; j < iterations_EM_; ++j)
            {
                // Likelihood of a datum given that it is an inlier
                for (size_t i = 0; i < indices_size; ++i)
                    p_inlier_prob[i] = gamma * exp (- (distances[i] * distances[i] ) / 2 * (sigma_ * sigma_) ) /
                        (sqrt (2 * M_PI) * sigma_);

                // Likelihood of a datum given that it is an outlier
                p_outlier_prob = (1 - gamma) / v;

                gamma = 0;
                for (size_t i = 0; i < indices_size; ++i)
                    gamma += p_inlier_prob [i] / (p_inlier_prob[i] + p_outlier_prob);
                gamma /= static_cast<double>(sac_model_->getIndices ()->size ());
            }

            // Find the log likelihood of the model -L = -sum [log (pInlierProb + pOutlierProb)]
            double d_cur_penalty = 0;
            for (size_t i = 0; i < indices_size; ++i)
                d_cur_penalty += log (p_inlier_prob[i] + p_outlier_prob);
            d_cur_penalty = - d_cur_penalty;

            // Better match ?
            if (d_cur_penalty < d_best_penalty)
            {
                d_best_penalty = d_cur_penalty;

                // Save the current model/coefficients selection as being the best so far
                model_              = selection;
                model_coefficients_ = model_coefficients;

                n_inliers_count = 0;
                // Need to compute the number of inliers for this model to adapt k
                for (size_t i = 0; i < distances.size (); ++i)
                    if (distances[i] <= 2 * sigma_)
                        n_inliers_count++;

                // Compute the k parameter (k=log(z)/log(1-w^n))
                double w = static_cast<double> (n_inliers_count) / static_cast<double> (sac_model_->getIndices ()->size ());
                double p_no_outliers = 1 - pow (w, static_cast<double> (selection.size ()));
                p_no_outliers = (std::max) (std::numeric_limits<double>::epsilon (), p_no_outliers);       // Avoid division by -Inf
                p_no_outliers = (std::min) (1 - std::numeric_limits<double>::epsilon (), p_no_outliers);   // Avoid division by 0.
                k = log (1 - probability_) / log (p_no_outliers);
            }

            ++iterations_;
            if (debug_verbosity_level > 1)
                SM_DEBUG ("[sm::MaximumLikelihoodSampleConsensus::computeModel] Trial %d out of %d. Best penalty is %f.\n", iterations_, static_cast<int> (ceil (k)), d_best_penalty);
            if (iterations_ > max_iterations_)
            {
                if (debug_verbosity_level > 0)
                    SM_DEBUG ("[sm::MaximumLikelihoodSampleConsensus::computeModel] MLESAC reached the maximum number of trials.\n");
                break;
            }
        }

        if (model_.empty ())
        {
            if (debug_verbosity_level > 0)
                SM_DEBUG ("[sm::MaximumLikelihoodSampleConsensus::computeModel] Unable to find a solution!\n");
            return (false);
        }

        // Iterate through the 3d points and calculate the distances from them to the model again
        sac_model_->getDistancesToModel (model_coefficients_, distances);
        std::vector<int> &indices = *sac_model_->getIndices ();
        if (distances.size () != indices.size ())
        {
            SM_ERROR ("[sm::MaximumLikelihoodSampleConsensus::computeModel] Estimated distances (%zu) differs than the normal of indices (%zu).\n", distances.size (), indices.size ());
            return (false);
        }

        inliers_.resize (distances.size ());
        // Get the inliers for the best model found
        n_inliers_count = 0;
        for (size_t i = 0; i < distances.size (); ++i)
            if (distances[i] <= 2 * sigma_)
                inliers_[n_inliers_count++] = indices[i];

        // Resize the inliers vector
        inliers_.resize (n_inliers_count);

        if (debug_verbosity_level > 0)
            SM_DEBUG ("[sm::MaximumLikelihoodSampleConsensus::computeModel] Model: %zu size, %d inliers.\n", model_.size (), n_inliers_count);

        return (true);
    }

    
} // namespace aslam
