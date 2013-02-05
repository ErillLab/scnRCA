import java.util.Comparator;


public class CDSComparatorLocation implements Comparator<CDS> {
	public CDSComparatorLocation()
	{
		
	}
	public int compare(CDS c1, CDS c2)
	{
		return c1.location() - c2.location();
	}
}
