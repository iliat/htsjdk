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

    public static void setDuplicateScore(final SAMRecord record, final ScoringStrategy scoringStrategy) {
        record.setAttribute(SAMTag.DS.name(), computeDuplicateScore(record, scoringStrategy));
    }

    public static void setDuplicateScore(final SAMRecord record1, final SAMRecord record2, final ScoringStrategy scoringStrategy) {
        final int duplicateScore = computeDuplicateScore(record1, record2, scoringStrategy);
        record1.setAttribute(SAMTag.DS.name(), duplicateScore);
        record2.setAttribute(SAMTag.DS.name(), duplicateScore);
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
    public static int computeDuplicateScore(final SAMRecord record, final ScoringStrategy scoringStrategy) {
        if (record.getReadPairedFlag()) {
            throw new SAMException("The read must be not be paired: " + record.getReadName());
        }
        return computeDuplicateScore(record, null, scoringStrategy);
    }

    /**
     * Returns the duplicate score computed from the given pair of records.  The pair of records are assumed to be mates.
     */
    public static int computeDuplicateScore(final SAMRecord record1, final SAMRecord record2, final ScoringStrategy scoringStrategy) {
        int score = 0;

        if (null != record2 && (!record1.getReadPairedFlag() || !record2.getReadPairedFlag())) {
            throw new SAMException("The read must be paired: " + record1.getReadName());
        }

        switch (scoringStrategy) {
            case SUM_OF_BASE_QUALITIES:
                score += getSumOfBaseQualities(record1);
                if (null != record2) score += getSumOfBaseQualities(record2);
                break;
            case TOTAL_MAPPED_REFERENCE_LENGTH:
                if (!record1.getReadUnmappedFlag()) {
                    score += record1.getCigar().getReferenceLength();
                }
                if (null != record2 && !record2.getReadUnmappedFlag()) {
                    score += record2.getCigar().getReferenceLength();
                }
                break;
        }
        return score;
    }

    /**
     * Compare two records based on their duplicate scores.  The duplicate scores for each record is assume to be
     * pre-computed by computeDuplicateScore and stored in the "DS" tag.  If the scores are equal, we break
     * ties based on mapping quality (added to the mate's mapping quality if paired and mapped), then library/read name.
     */
    /** We allow different scoring strategies. We return <0 if rec1 has a better strategy than rec2. */
    public static int compare(final SAMRecord rec1, final SAMRecord rec2) {
        int cmp;

        // always prefer paired over non-paired
        if (rec1.getReadPairedFlag() != rec2.getReadPairedFlag()) return rec1.getReadPairedFlag() ? 1 : -1;

        // Get the primary duplicate score
        Integer duplicateScore1 = (Integer)rec1.getAttribute(SAMTagUtil.getSingleton().DS);
        if (null == duplicateScore1) throw new SAMException("DS tag not found for: " + rec1.getReadName());
        Integer duplicateScore2 = (Integer)rec2.getAttribute(SAMTagUtil.getSingleton().DS);
        if (null == duplicateScore2) throw new SAMException("DS tag not found for: " + rec2.getReadName());
        cmp = duplicateScore1 - duplicateScore2;

        // Get the mapping quality
        // NB: must have the mate's mapping quality too!
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

        /**
         * Finally, use library ID and read name
         * This is important because we cannot control the order in which reads appear for reads that are comparable up to now (i.e. cmp == 0).  We want to deterministically
         * choose them, and so we need this.
        */
        if (0 == cmp) cmp = SAMUtils.getCanonicalRecordName(rec2).compareTo(SAMUtils.getCanonicalRecordName(rec1)); // since we will return -cmp, reverse the comparison

        return -cmp;
    }

}
