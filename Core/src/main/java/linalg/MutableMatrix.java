package linalg;

public interface MutableMatrix extends Matrix, MutableTensor
{
	void deleteRow(int row);
	void deleteColumn(int column);
	
}
