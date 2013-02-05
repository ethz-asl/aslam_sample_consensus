
namespace aslam {

    template<typename P>
    SampleConsensus<P>::SampleConsensus(int maxIterations, double threshold, double probability) :
        max_iterations_(maxIterations),
        threshold_(threshold),
        probability_(probability)
    {

    }

    template<typename P>
    SampleConsensus<P>::~SampleConsensus(){}

}
