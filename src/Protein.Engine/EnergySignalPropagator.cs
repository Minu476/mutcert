using System;

namespace Protein.Engine;

/// <summary>
/// Applies the monotone-narrowing interval update rule defined in the MutCert spec (v3.3):
///
///   lo_new = max(lo_prev, DDG_k − ε₀ / √k)
///   hi_new = min(hi_prev, DDG_k + ε₀ / √k)
///
/// where DDG_k is the running cumulative sum of energy deltas after k BFS steps,
/// and ε₀ is a per-family calibration parameter (phase 1: training-set standard deviation).
///
/// Guarantees: lo is non-decreasing and hi is non-increasing across steps unless a
/// crossover (lo > hi) is detected, in which case the interval resets to the
/// unconstrained current estimate and CrossoverOccurred is flagged for one step.
/// </summary>
public class EnergySignalPropagator
{
    private readonly double _epsilon0;
    private double _cumulativeDDG;
    private double _lo;
    private double _hi;
    private readonly bool _hasInitial;

    public double Lo            => _lo;
    public double Hi            => _hi;
    public double Width         => _hi - _lo;
    public double CumulativeDDG => _cumulativeDDG;
    public int    StepCount     { get; private set; }

    /// <summary>
    /// True if the most recent Apply() call triggered a crossover recovery.
    /// Resets to false on the next Apply() call.
    /// </summary>
    public bool CrossoverOccurred { get; private set; }

    /// <summary>
    /// Default constructor. The interval is set by the first Apply() call:
    /// [DDG_1 − ε₀, DDG_1 + ε₀].
    /// </summary>
    public EnergySignalPropagator(double epsilon0)
    {
        if (epsilon0 <= 0)
            throw new ArgumentOutOfRangeException(nameof(epsilon0), "ε₀ must be positive.");

        _epsilon0       = epsilon0;
        _cumulativeDDG  = 0.0;
        _lo             = double.NegativeInfinity;
        _hi             = double.PositiveInfinity;
        _hasInitial     = false;
        StepCount       = 0;
        CrossoverOccurred = false;
    }

    /// <summary>
    /// Spec-correct constructor (Phase 6).  The initial interval is pre-seeded from the
    /// raw Miyazawa-Jernigan site estimate computed before BFS begins:
    ///   [lo₀, hi₀] = [ΔΔG_MJ − ε₀, ΔΔG_MJ + ε₀]
    ///
    /// The running <see cref="CumulativeDDG"/> also starts at <paramref name="initialDdgMj"/>
    /// so subsequent BFS deltas accumulate on top of the pre-computed site estimate.
    /// </summary>
    public EnergySignalPropagator(double epsilon0, double initialDdgMj)
    {
        if (epsilon0 <= 0)
            throw new ArgumentOutOfRangeException(nameof(epsilon0), "ε₀ must be positive.");

        _epsilon0       = epsilon0;
        _cumulativeDDG  = initialDdgMj;
        _lo             = initialDdgMj - epsilon0;
        _hi             = initialDdgMj + epsilon0;
        _hasInitial     = true;
        StepCount       = 0;
        CrossoverOccurred = false;
    }

    /// <summary>
    /// Accepts one BFS energy delta and narrows the [lo, hi] interval.
    /// Call once per BFS propagation step, in step order.
    /// </summary>
    public void Apply(double delta)
    {
        CrossoverOccurred = false;

        StepCount++;
        _cumulativeDDG += delta;

        double uncertainty  = _epsilon0 / Math.Sqrt(StepCount);
        double proposedLo   = _cumulativeDDG - uncertainty;
        double proposedHi   = _cumulativeDDG + uncertainty;

        // When no pre-seed is given, step 1 sets the interval directly.
        // When a pre-seed (initialDdgMj) is given, max/min narrows from step 1.
        if (StepCount == 1 && !_hasInitial)
        {
            _lo = proposedLo;
            _hi = proposedHi;
        }
        else
        {
            _lo = Math.Max(_lo, proposedLo);
            _hi = Math.Min(_hi, proposedHi);
        }

        // Crossover recovery: conflicting BFS steps can invert the interval.
        // Reset to the unconstrained current estimate so the run can continue.
        if (_lo > _hi)
        {
            CrossoverOccurred = true;
            _lo = proposedLo;
            _hi = proposedHi;
        }
    }
}
