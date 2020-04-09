# Algorithms

```@meta
CurrentModule = BioSimulator
```

## Exact algorithms

### Direct methods

```@docs
Direct
EnhancedDirect
SortingDirect
```

### First reaction methods

```@docs
FirstReaction
NextReaction
```

### Rejection methods

```@docs
RejectionSSA
```

!!! warning

    This implementation needs additional testing and is not guaranteed to be optimized.

## Approximate algorithms

!!! warning

    These implementations are a work in progress.
    The `TauLeapingDG2001` and `TauLeapingDGLP2003` methods need "hybrid" counterparts.

```@docs
TauLeapingDG2001
TauLeapingDGLP2003
StepAnticipation
HybridSAL
```
