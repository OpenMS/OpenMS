//! [MRMTransitionGroup]

typedef MRMTransitionGroup<MSChromatogram, ReactionMonitoringTransition> TrGroup;
  TrGroup createTransitionGroup()
  {
    TrGroup tr_group;
    tr_group.addChromatogram(MSChromatogram(), “transition1”);
    tr_group.addTransition(ReactionMonitoringTransition(), “transition1”);
    tr_group.addChromatogram(MSChromatogram(), “transition2”);
    tr_group.addTransition(ReactionMonitoringTransition(), “transition2”);
    tr_group.setTransitionGroupID(“tr_peptideA”);
    return tr_group;
  }

//! [MRMTransitionGroup]
