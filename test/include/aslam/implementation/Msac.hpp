
namespace aslam {

    template<typename P>
    Msac<P>::Msac(int maxIterations, double threshold, double probability) :
        SampleConsensus<P>(maxIterations, threshold, probability)
    {

    }

    template<typename P>
    Msac<P>::~Msac(){}


    template<typename PROBLEM_T>
    bool Msac<PROBLEM_T>::computeModel(int debug_verbosity_level)
    {
        typedef PROBLEM_T problem_t;
        typedef typename problem_t::model_t model_t;

        // Warn and exit if no threshold was set
        if (threshold_ == std::numeric_limits<double>::max())
        {
            SM_ERROR ("[sm::MEstimatorSampleConsensus::computeModel] No threshold set!\n");
            return (false);
        }

        iterations_ = 0;
        double d_best_penalty = std::numeric_limits<double>::max();
        double k = 1.0;

        std::vector<int> best_model;
        std::vector<int> selection;
        model_t model_coefficients;
        std::vector<double> distances;

        int n_inliers_count = 0;
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

            double d_cur_penalty = 0;
            // Iterate through the 3d points and calculate the distances from them to the model
            sac_model_->getDistancesToModel (model_coefficients, distances);
    
            if (distances.empty () && k > 1.0)
                continue;

            for (size_t i = 0; i < distances.size (); ++i)
                d_cur_penalty += (std::min) (distances[i], threshold_);

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
                    if (distances[i] <= threshold_)
                        ++n_inliers_count;

                // Compute the k parameter (k=log(z)/log(1-w^n))
                double w = static_cast<double> (n_inliers_count) / static_cast<double> (sac_model_->getIndices ()->size ());
                double p_no_outliers = 1.0 - pow (w, static_cast<double> (selection.size ()));
                p_no_outliers = (std::max) (std::numeric_limits<double>::epsilon (), p_no_outliers);       // Avoid division by -Inf
                p_no_outliers = (std::min) (1.0 - std::numeric_limits<double>::epsilon (), p_no_outliers);   // Avoid division by 0.
                k = log (1.0 - probability_) / log (p_no_outliers);
            }

            ++iterations_;
            if (debug_verbosity_level > 1)
                SM_DEBUG ("[sm::MEstimatorSampleConsensus::computeModel] Trial %d out of %d. Best penalty is %f.\n", iterations_, static_cast<int> (ceil (k)), d_best_penalty);
            if (iterations_ > max_iterations_)
            {
                if (debug_verbosity_level > 0)
                    SM_DEBUG ("[sm::MEstimatorSampleConsensus::computeModel] MSAC reached the maximum number of trials.\n");
                break;
            }
        }

        if (model_.empty ())
        {
            if (debug_verbosity_level > 0)
                SM_DEBUG ("[sm::MEstimatorSampleConsensus::computeModel] Unable to find a solution!\n");
            return (false);
        }

        // Iterate through the 3d points and calculate the distances from them to the model again
        sac_model_->getDistancesToModel (model_coefficients_, distances);
        std::vector<int> &indices = *sac_model_->getIndices ();

        if (distances.size () != indices.size ())
        {
            SM_ERROR ("[sm::MEstimatorSampleConsensus::computeModel] Estimated distances (%zu) differs than the normal of indices (%zu).\n", distances.size (), indices.size ());
            return (false);
        }

        inliers_.resize (distances.size ());
        // Get the inliers for the best model found
        n_inliers_count = 0;
        for (size_t i = 0; i < distances.size (); ++i)
            if (distances[i] <= threshold_)
                inliers_[n_inliers_count++] = indices[i];

        // Resize the inliers vector
        inliers_.resize (n_inliers_count);

        if (debug_verbosity_level > 0)
            SM_DEBUG ("[sm::MEstimatorSampleConsensus::computeModel] Model: %zu size, %d inliers.\n", model_.size (), n_inliers_count);

        return (true);

    }

    
} // namespace aslam
