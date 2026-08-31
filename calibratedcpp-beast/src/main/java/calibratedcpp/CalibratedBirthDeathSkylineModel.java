package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.Tensor;

import java.util.*;

/**
 * @author Marcus Overwater
 */

@Description("Extension of CalibratedCoalescentPointProcess, implements node age density and CDF of " +
        "for the BDSKY model with piecewise constant birth and death rates with incomplete extant sampling.")
public class CalibratedBirthDeathSkylineModel extends CalibratedCoalescentPointProcess {
    public Input<SkylineParameter> birthRateInput =
            new Input<>("birthRate", "Skyline parameter for the birthRate (λ)", (SkylineParameter) null);
    public Input<SkylineParameter> deathRateInput =
            new Input<>("deathRate", "Skyline parameter for the deathRate (µ)", (SkylineParameter) null);
    public Input<SkylineParameter> diversificationRateInput =
            new Input<>("diversificationRate", "Skyline parameter for the diversificationRate (λ - µ)", (SkylineParameter) null);
    public Input<SkylineParameter> reproductiveNumberInput =
            new Input<>("reproductiveNumber", "Skyline parameter for the reproductiveNumber (λ/µ)", (SkylineParameter) null);
    public Input<SkylineParameter> turnoverInput =
            new Input<>("turnover", "Skyline parameter for the turnover (µ/λ)", (SkylineParameter) null);
    public Input<RealScalar<UnitInterval>> samplingProbabilityInput =
            new Input<>("rho", "Sampling probability (⍴)", (RealScalar<UnitInterval>) null);


    protected double[] intervalStartTimes, lambda, r, cumulativeIntegral, cumulativeExpR;
    protected double rho;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // 1. Store the PARENT inputs in the array, not the children yet.
        Input<SkylineParameter>[] skylineInputs = new Input[]{
                birthRateInput, deathRateInput, reproductiveNumberInput,
                diversificationRateInput, turnoverInput
        };

        int specifiedRates = 0;
        StringBuilder whichSpecified = new StringBuilder();

        for (int i = 0; i < 5; i++) {
            SkylineParameter sp = skylineInputs[i].get();

            // Check if the SkylineParameter itself exists
            if (sp != null) {
                // Now it is safe to access fields
                Tensor<?, ?> rateP = sp.valuesInput.get();
                Tensor<?, ?> timeP = sp.changeTimesInput.get();

                  if (rateP != null) {
                    specifiedRates++;
                    if (!whichSpecified.isEmpty()) whichSpecified.append(", ");
                    whichSpecified.append(skylineInputs[i].getName());

                    if (timeP != null) {
                        int rateDim = rateP.size();
                        int timeDim = timeP.size();
                        if (rateDim != timeDim + 1) {
                            throw new IllegalArgumentException("Dimension mismatch for " + skylineInputs[i].getName() +
                                    ": Explicit change times provided (" + timeDim + "), so rate dimension must be (" + (timeDim + 1) + "). Found: " + rateDim);
                        }
                    }
                }
            }
        }

        if (samplingProbabilityInput.get() == null) {
            throw new IllegalArgumentException("rho parameter must be specified.");
        }

        if (reproductiveNumberInput.get() != null && turnoverInput.get() != null)
            throw new IllegalArgumentException("Cannot specify both reproductiveNumber and turnover.");

        if (samplingProbabilityInput.get() == null)
            throw new IllegalArgumentException("Sampling probability (rho) must be specified.");

        if (specifiedRates != 2)
            throw new IllegalArgumentException("Exactly TWO rates must be specified. Found " + specifiedRates + " (" + whichSpecified + ")");

        updateIntervals(true);
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        if (!updateIntervals(false)) {
            return Double.NEGATIVE_INFINITY;
        }
        return super.calculateTreeLogLikelihood(tree);
    }

    private boolean updateIntervals(boolean strict) {
        super.updateModel();

        rho = samplingProbabilityInput.get().get();
        double rootHeight = treeInput.get().getRoot().getHeight();
        double originVal = !conditionOnRoot ? originInput.get().get() : 0.0;
        double maxT = maxTime;

        if (!conditionOnRoot && rootHeight >= originVal) return false;

        // 1. Process Times safely
        List<Double> bTimes = processSafe(birthRateInput, maxT);
        List<Double> dTimes = processSafe(deathRateInput, maxT);
        List<Double> rTimes = processSafe(reproductiveNumberInput, maxT);
        List<Double> divTimes = processSafe(diversificationRateInput, maxT);
        List<Double> tTimes = processSafe(turnoverInput, maxT);

        // 2. Create Master Timeline (Union of all times).
        // Only cuts falling strictly inside (0, maxT) become integration boundaries: a cut outside
        // that range bounds an epoch that lies wholly outside the process and so contributes
        // nothing. The per-parameter lists above deliberately retain those out-of-range cuts
        SortedSet<Double> timesSet = new TreeSet<>();
        timesSet.add(0.0);
        for (List<Double> cuts : List.of(bTimes, dTimes, rTimes, divTimes, tTimes))
            for (double t : cuts)
                if (t > 1e-10 && t < maxT - 1e-10) timesSet.add(t);

        intervalStartTimes = timesSet.stream().mapToDouble(d -> d).toArray();
        int n = intervalStartTimes.length;
        lambda = new double[n];
        r = new double[n];
        cumulativeIntegral = new double[n];
        cumulativeExpR = new double[n];

        double logRunningSum = Double.NEGATIVE_INFINITY;
        double rRunningSum = 0.0;
        cumulativeIntegral[0] = Double.NEGATIVE_INFINITY;
        cumulativeExpR[0] = 0.0;

        for (int j = 0; j < n; j++) {
            double t = intervalStartTimes[j];

            // 2. Get Values Safely
            double vL = getValSafe(birthRateInput, bTimes, t);
            double vM = getValSafe(deathRateInput, dTimes, t);
            double vR = getValSafe(reproductiveNumberInput, rTimes, t);
            double vD = getValSafe(diversificationRateInput, divTimes, t);
            double vT = getValSafe(turnoverInput, tTimes, t);

            double l = 0, m = 0;

            if (birthRateInput.get() != null && deathRateInput.get() != null) { l = vL; m = vM; }
            else if (birthRateInput.get() != null && diversificationRateInput.get() != null) { l = vL; m = vL - vD; }
            else if (birthRateInput.get() != null && reproductiveNumberInput.get() != null) { l = vL; m = vL / vR; }
            else if (birthRateInput.get() != null && turnoverInput.get() != null) { l = vL; m = vL * vT; }
            else if (deathRateInput.get() != null && diversificationRateInput.get() != null) { m = vM; l = vM + vD; }
            else if (deathRateInput.get() != null && turnoverInput.get() != null) { m = vM; l = vM / vT; }
            else if (deathRateInput.get() != null && reproductiveNumberInput.get() != null) { m = vM; l = m * vR;}
            else if (diversificationRateInput.get() != null && reproductiveNumberInput.get() != null) { m = vD / (vR - 1.0); l = m * vR; }
            else if (diversificationRateInput.get() != null && turnoverInput.get() != null) { l = vD / (1.0 - vT); m = l * vT; }

            if (l < 0.0 || m < 0.0) {
                if (strict) throw new IllegalArgumentException(String.format(
                        "%s: the supplied parameterisation (%s) yields a negative rate in interval %d "
                                + "(age %.6g): birth λ = %.6g, death µ = %.6g. Adjust the starting values or priors "
                                + "so both birth and death rates stay non-negative.",
                        getID(), activeParameterisation(), j, t, l, m));
                return false;
            }

            lambda[j] = l;
            r[j] = l - m;

            if (j < n - 1) {
                double dt = intervalStartTimes[j+1] - t;
                logRunningSum = logSumExp(logRunningSum, calculateSegment(lambda[j], r[j], dt, rRunningSum));
                rRunningSum += (r[j] * dt);
                cumulativeIntegral[j+1] = logRunningSum;
                cumulativeExpR[j+1] = rRunningSum;
            }
        }
        return true;
    }

    /**
     * Processes input times, converting them all to "Age" (Time from Present).
     * @param reverse If TRUE: Input is already Age (Distance from Present).
     * If FALSE: Input is Distance from Root (needs conversion).
     */
    private List<Double> processInput(Tensor<?, ?> rateP, Tensor<?, ?> timeP,
                                      boolean relative, boolean reverse, double maxTime) {
        List<Double> times = new ArrayList<>();
        if (rateP == null) return times;

        if (timeP != null) {
            // --- Explicit Change Times ---
            // No sorting: SkylineParameter.initAndValidate() requires changeTimes to be strictly
            // increasing, and each conversion below is monotone, so the resulting ages are already
            // ordered — ascending when reverse=TRUE, descending when reverse=FALSE (both from-root
            // conversions negate v), which the single reverse() below corrects. Sorting instead
            // would silently reinterpret an out-of-order parameter rather than letting it error.
            double[] rawTimes = new double[timeP.size()];
            for (int i = 0; i < timeP.size(); i++) rawTimes[i] = (Double) timeP.get(i);
            for (double v : rawTimes) {
                double tAge;
                if (relative) {
                    // relative [0,1]
                    // If reverse=TRUE (Age): v=0 is Present. Age = v * maxTime.
                    // If reverse=FALSE (FromRoot): v=0 is Root. Age = maxTime * (1-v).
                    tAge = reverse ? (v * maxTime) : (maxTime * (1.0 - v));
                } else {
                    // Absolute
                    // If reverse=TRUE (Age): Input v is Age.
                    // If reverse=FALSE (FromRoot): Input v is time since root. Age = maxTime - v.
                    tAge = reverse ? v : (maxTime - v);
                }

                // Out-of-range cuts are deliberately retained, so that this list always has
                // exactly rateP.size()-1 entries and getVal()'s index mapping stays exact.
                // Change times are validated non-negative and strictly increasing (and <= 1 when
                // relative), so the only way out of the process is an absolute time past maxTime —
                // always the trailing k. Since values share the direction of the times, that
                // strands values[size-k .. size-1]: those epochs lie wholly outside [0, maxTime]
                // and correctly contribute nothing. (Under !timesAreAges those k cuts convert to
                // negative ages and land at the front of this list after the reverse() below.)
                times.add(tAge);
            }
            if (!reverse) Collections.reverse(times);
        } else {
            // --- Implicit Equidistant Times ---
            // If explicit times are missing, we assume equidistant intervals over maxTime.
            int numIntervals = rateP.size();
            if (numIntervals > 1) {
                double width = maxTime / numIntervals;
                for (int i = 1; i < numIntervals; i++) {
                    times.add(width * i);
                }
            }
        }
        return times;
    }

    /**
     * Retrieves the rate value for a given time t.
     * * @param t The current time (Age).
     */
    private double getVal(Tensor<?, ?> p, List<Double> cuts, double t, boolean reverse) {
        if (p == null) return 0;

        // Find which time interval 't' falls into based on the cuts (which are Ages).
        // If t is in [0, t1], index is 0.
        int idx = Collections.binarySearch(cuts, t);
        int intervalIndex = (idx >= 0) ? idx + 1 : -idx - 1;

        // ORIENTATION LOGIC: values[0] is always in [0,changeTime[0]]
        // If timesAreAges, values are specified Present-to-Root:
        //   values[0]    in [0, cut[0]]            (interval before the first change time)
        //   values[Last] in [cut[Last], maxTime]
        // If !timesAreAges, values are specified Root-to-Present, and the mapping is mirrored:
        //   values[0]    in [cut[Last], maxTime]
        //   values[Last] in [0, cut[0]]
        int pIdx = reverse ? intervalIndex : p.size() - 1 - intervalIndex;

        return (Double) p.get(Math.max(0, Math.min(pIdx, p.size() - 1)));
    }

    private String activeParameterisation() {
        List<String> set = new ArrayList<>();
        if (birthRateInput.get()          != null) set.add("birthRate");
        if (deathRateInput.get()          != null) set.add("deathRate");
        if (diversificationRateInput.get()!= null) set.add("diversificationRate");
        if (reproductiveNumberInput.get() != null) set.add("reproductiveNumber");
        if (turnoverInput.get()           != null) set.add("turnover");
        return String.join(" + ", set);
    }

    // Helper to extract times without crashing on null inputs
    private List<Double> processSafe(Input<SkylineParameter> input, double maxT) {
        SkylineParameter sp = input.get();
        if (sp == null) return new ArrayList<>();
        return processInput(sp.valuesInput.get(), sp.changeTimesInput.get(),
                sp.isRelative, sp.isReverse, maxT);
    }

    // Helper to get value without crashing on null inputs
    private double getValSafe(Input<SkylineParameter> input, List<Double> times, double t) {
        SkylineParameter sp = input.get();
        if (sp == null || sp.valuesInput.get() == null) return 0.0;
        return getVal(sp.valuesInput.get(), times, t, sp.timesAreAgesInput.get());
    }

    // --- Math Helpers ---

    private double calculateSegment(double l, double r_val, double dt, double expOffset) {
        double logTerm;
        double x = r_val * dt;
        if (Math.abs(r_val) < 1e-9) {
            logTerm = Math.log(rho) + Math.log(l) + Math.log(dt);
        } else if (r_val > 0) {
            double logDiff = Math.log(Math.expm1(x));
            logTerm = Math.log(rho) + Math.log(l) - Math.log(r_val) + logDiff;
        } else {
            double logDiff = Math.log(-Math.expm1(x));
            logTerm = Math.log(rho) + Math.log(l) - Math.log(-r_val) + logDiff;
        }
        return logTerm + expOffset;
    }

    @Override
    public double calculateLogNodeAgeDensity(double time) {
        int m = getInterval(time);
        double dtP = time - intervalStartTimes[m];
        double logY = logSumExp(0.0, logSumExp(cumulativeIntegral[m], calculateSegment(lambda[m], r[m], dtP, cumulativeExpR[m])));
        return (Math.log(rho) + Math.log(lambda[m]) + cumulativeExpR[m] + (r[m] * dtP)) - 2.0 * logY;
    }

    @Override
    public double calculateLogNodeAgeCDF(double time) {
        int m = getInterval(time);
        double logInt = logSumExp(cumulativeIntegral[m], calculateSegment(lambda[m], r[m], time - intervalStartTimes[m], cumulativeExpR[m]));
        return logInt - logSumExp(0.0, logInt);
    }

    /**
     * Computes log(1 - Q(time)) directly rather than via the CDF.
     */
    @Override
    public double calculateLogNodeAgeSurvival(double time) {
        int m = getInterval(time);
        double logInt = logSumExp(cumulativeIntegral[m], calculateSegment(lambda[m], r[m], time - intervalStartTimes[m], cumulativeExpR[m]));
        return -logSumExp(0.0, logInt);
    }

    private int getInterval(double t) {
        int i = Arrays.binarySearch(intervalStartTimes, t);
        return i < 0 ? Math.max(0, -i - 2) : i;
    }

    private double logSumExp(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b; if (b == Double.NEGATIVE_INFINITY) return a;
        return Math.max(a, b) + Math.log1p(Math.exp(-Math.abs(a - b)));
    }
}