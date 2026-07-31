import calibratedcpp.lphy.tree.CPPTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * The uncalibrated tip count is determined, in priority order, by explicit n, then the number of
 * supplied taxa, then random. So supplying taxa fixes both the count and the names, and random-N
 * only occurs when no taxa are given (making the taxa-vs-random-count mismatch impossible).
 */
public class TaxaImpliesNTest {

    private static void collectLeafIds(TimeTreeNode node, Set<String> ids) {
        List<TimeTreeNode> children = node.getChildren();
        if (children == null || children.isEmpty()) { ids.add(node.getId()); return; }
        for (TimeTreeNode child : children) collectLeafIds(child, ids);
    }

    private CPPTree cpp(String[] taxa, Integer n, Double rootAge) {
        return new CPPTree(new Value<>("", 1.5), new Value<>("", 0.5), null, null, new Value<>("", 1.0),
                taxa == null ? null : new Value<>("", taxa),
                n == null ? null : new Value<>("", n),
                rootAge == null ? null : new Value<>("", rootAge), null);
    }

    @Test
    public void taxaDetermineTipCountAndNames() {
        Set<String> expected = Set.of("A", "B", "C", "D");

        // taxa, no n, fixed root -> 4 tips named A..D
        Set<String> ids = new HashSet<>();
        collectLeafIds(cpp(new String[]{"A", "B", "C", "D"}, null, 5.0).sample().value().getRoot(), ids);
        assertEquals(expected, ids, "taxa determine tip names when n is unset (fixed root)");

        // taxa, no n, no rootAge -> random stem, still 4 tips named A..D
        Set<String> ids2 = new HashSet<>();
        collectLeafIds(cpp(new String[]{"A", "B", "C", "D"}, null, null).sample().value().getRoot(), ids2);
        assertEquals(expected, ids2, "taxa determine tips in the random-stem case too");
    }

    @Test
    public void taxaLongerThanExplicitNThrows() {
        CPPTree g = cpp(new String[]{"A", "B", "C", "D"}, 3, 5.0);
        assertThrows(IllegalArgumentException.class, g::sample);
    }
}
