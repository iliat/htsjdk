package htsjdk.samtools;

/**
 * This class helps us compute and compare duplicate scores, which are used for selecting the non-duplicate
 * during duplicate marking (see MarkDuplicates).
 * @author nhomer
 */
public class DuplicateScoringStrategy {

    public enum ScoringStrategy {
        SUM_OF_BASE_QUALITIES,
        TOTAL_MAPPED_REFERENCE_LENGTH
    }

    /** Calculates a score for the read which is the sum of scores over Q15. */
    private static short getSumOfBaseQualities(final SAMRecord rec) {
        short score = 0;
        for (final byte b : rec.getBaseQualities()) {
            if (b >= 15) score += b;
        }

        return score;
    }

    /**
     * Returns the duplicate score computed from the given fragment.
     */
    public static short computeDuplicateScore(final SAMRecord record, final ScoringStrategy scoringStrategy) {
        short score = 0;

        switch (scoringStrategy) {
            case SUM_OF_BASE_QUALITIES:
                score += getSumOfBaseQualities(record);
                break;
            case TOTAL_MAPPED_REFERENCE_LENGTH:
                if (!record.getReadUnmappedFlag()) {
                    score += record.getCigar().getReferenceLength();
                }
                break;
        }
        return score;
    }

    /**
     * Returns the duplicate score computed from the given fragment.
     *
     * If true is given to assumeMateCigar, then any score that can use the mate cigar to to compute the mate's score will return the score
     * computed on both ends.
     */
    public static short computeDuplicateScore(final SAMRecord record, final ScoringStrategy scoringStrategy, final boolean assumeMateCigar) {
        short score = 0;

        switch (scoringStrategy) {
            case SUM_OF_BASE_QUALITIES:
                score += getSumOfBaseQualities(record);
                break;
            case TOTAL_MAPPED_REFERENCE_LENGTH:
                if (!record.getReadUnmappedFlag()) {
                    score += record.getCigar().getReferenceLength();
                }
                if (assumeMateCigar && record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {
                    score += SAMUtils.getMateCigar(record).getReferenceLength();
                }
                break;
        }
        return score;
    }

    /**
     * Compare two records based on their duplicate scores.  The duplicate scores for each record is assume to be
     * pre-computed by computeDuplicateScore and stored in the "DS" tag.  If the scores are equal, we break
     * ties based on mapping quality (added to the mate's mapping quality if paired and mapped), then library/read name.
     *
     * If true is given to assumeMateCigar, then any score that can use the mate cigar to to compute the mate's score will return the score
     * computed on both ends.
     *
     * We allow different scoring strategies. We return <0 if rec1 has a better strategy than rec2.
     */
    public static int compare(final SAMRecord rec1, final SAMRecord rec2, final ScoringStrategy scoringStrategy, final boolean assumeMateCigar) {
        int cmp;

        // always prefer paired over non-paired
        if (rec1.getReadPairedFlag() != rec2.getReadPairedFlag()) return rec1.getReadPairedFlag() ? 1 : -1;

        // Get the primary duplicate score
        // TODO: remove me
        /*
        Integer duplicateScore1 = (Integer)rec1.getAttribute(SAMTagUtil.getSingleton().DS);
        if (null == duplicateScore1) throw new SAMException("DS tag not found for: " + rec1.getReadName());
        Integer duplicateScore2 = (Integer)rec2.getAttribute(SAMTagUtil.getSingleton().DS);
        if (null == duplicateScore2) throw new SAMException("DS tag not found for: " + rec2.getReadName());
        cmp = duplicateScore1 - duplicateScore2;
        */
        cmp = computeDuplicateScore(rec1, scoringStrategy, assumeMateCigar) - computeDuplicateScore(rec2, scoringStrategy, assumeMateCigar);

        // Get the mapping quality
        // NB: must have the mate's mapping quality too!
        /*
        if (0 == cmp) {
            Integer mateMappingQuality1 = 0;
            if (rec1.getReadPairedFlag() && !rec1.getMateUnmappedFlag()) {
                mateMappingQuality1 = (Integer)rec1.getAttribute(SAMTagUtil.getSingleton().MQ);
                if (null == mateMappingQuality1) throw new SAMException("Missing MQ for read: " + rec1.getReadName());
            }
            Integer mateMappingQuality2 = 0;
            if (rec2.getReadPairedFlag() && !rec2.getMateUnmappedFlag()) {
                mateMappingQuality2 = (Integer)rec2.getAttribute(SAMTagUtil.getSingleton().MQ);
                if (null == mateMappingQuality2) throw new SAMException("Missing MQ for read: " + rec2.getReadName());
            }
            // TODO: should we check that rec1/rec2 are themselves mapped?
            cmp = rec1.getMappingQuality() - rec2.getMappingQuality() + mateMappingQuality1 - mateMappingQuality2;
        }
        */

        /**
         * Finally, use library ID and read name
         * This is important because we cannot control the order in which reads appear for reads that are comparable up to now (i.e. cmp == 0).  We want to deterministically
         * choose them, and so we need this.
         */
        if (0 == cmp) cmp = SAMUtils.getCanonicalRecordName(rec2).compareTo(SAMUtils.getCanonicalRecordName(rec1)); // since we will return -cmp, reverse the comparison

        return -cmp;
    }

    /**
     * Compare two records based on their duplicate scores.  The duplicate scores for each record is assume to be
     * pre-computed by computeDuplicateScore and stored in the "DS" tag.  If the scores are equal, we break
     * ties based on mapping quality (added to the mate's mapping quality if paired and mapped), then library/read name.
     *
     * We allow different scoring strategies. We return <0 if rec1 has a better strategy than rec2.
     */
    public static int compare(final SAMRecord rec1, final SAMRecord rec2, final ScoringStrategy scoringStrategy) {
        return compare(rec1, rec2, scoringStrategy, false);
    }

}
