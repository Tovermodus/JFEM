package basic;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class StatLogger
{
	public static void log(final String message)
	{
		
		final BufferedWriter writer;
		try
		{
			writer = new BufferedWriter(new FileWriter("../dlm/statlog", true));
			writer.write(message + "\n");
			writer.close();
		} catch (final IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public static void clear()
	{
		final File f = new File("../dlm/statlog");
		if (f.exists())
			f.delete();
	}
	
	public static void log(final List<String> messages)
	{
		
		final BufferedWriter writer;
		try
		{
			writer = new BufferedWriter(new FileWriter("../dlm/statlog", true));
			for (final var message : messages)
				writer.write(message + "\n");
			writer.close();
		} catch (final IOException e)
		{
			e.printStackTrace();
		}
	}
}
