package basic;

import javax.swing.*;
import java.awt.*;

public class Interruptor
	implements Runnable
{
	private JFrame f;
	public volatile boolean running = true;
	
	public boolean isRunning()
	{
		return running;
	}
	
	@Override
	public void run()
	{
		f = new JFrame("Interrupt Solver");
		final JButton b = new JButton("Interrupt!");
		f.setLayout(new GridLayout());
		f.add(b);
		f.setBounds(0, 0, 300, 300);
		f.setVisible(true);
		b.addActionListener(e ->
		                    {
			                    System.out.println("Interrupt!!!");
			                    running = false;
			                    f.setVisible(false);
			                    f.dispose();
		                    });
		while (running)
		{
			try
			{
				Thread.sleep(30);
			} catch (final InterruptedException e)
			{
				e.printStackTrace();
			}
		}
		f.dispose();
	}
}
