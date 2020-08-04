package mixed;
import basic.Cell;
import basic.Face;

public interface MixedFESpace<CT extends Cell<CT,FT,ST>, FT extends Face<CT,FT,ST>,ST extends MixedShapeFunction<CT,
	FT,ST>>
{

}
