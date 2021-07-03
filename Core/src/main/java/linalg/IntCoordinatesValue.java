package linalg;

public class IntCoordinatesValue
{
	public final IntCoordinates coordinates;
	public final double value;
	public final boolean add;
	
	public IntCoordinatesValue(IntCoordinates coordinates, double value, boolean add)
	{
		this.coordinates = coordinates;
		this.value = value;
		this.add = add;
	}
}