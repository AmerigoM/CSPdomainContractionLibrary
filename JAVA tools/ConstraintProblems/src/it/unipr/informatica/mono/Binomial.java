package it.unipr.informatica.mono;

/*
	@file Binomial.java
	@author Amerigo Mancino, amerigo.mancino@gmail.com
	@version 1.0
*/

/** Implementazione del calcolo del coefficiente
  binomiale matematico */
public class Binomial {
	
    public static long binomial(int n, int k) {
    	if (k > n) {
    		System.out.println("Errore!");
    	    return -1;
    	}
        if (k > n-k)
            k = n-k;
 
        long b = 1;
        for (int i = 1, m = n; i <= k; i++, m--)
            b = b*m/i;
        
        return b;
    }

}
