# masters-thesis
In this thesis, the goal was to create a system on hardware capable of detecting a brain injury. For that purpose four features which can reveal an injury and are called Relative Power, Amplitude Assymetry, Coherence and Phase Difference, were implemented. The combination of the four features in a discriminant model i.e. a Discriminant Function(DF) is the tool that will reveal the presence or not of an injury. First, a MATLAB model was built for the features and the DF, which was then developed on hardware. The four feautures were realised as HLS IP accelerators and the DF was realised as a user level application developed with the SDK tool. Each accelerator implements the AXI4-stream protocol and is fed with data using DMAs. Also, the Coherence and Phase Difference features were combined on a signle accelerator on hardware in order to reduce the number of DMA instances that were used for this application. The target board was a ZCU104.
