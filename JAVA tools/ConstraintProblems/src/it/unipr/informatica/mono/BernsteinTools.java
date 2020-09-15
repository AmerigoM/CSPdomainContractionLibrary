package it.unipr.informatica.mono;

import Jama.*;

/** Implementazione di metodi per la conversione da
 * forma polinomiale a base di Bernstein */
public class BernsteinTools {

	/**********************************************************/
	/** Metodo per la conversione da forma di potenze in forma
	 * di Bernstein per polinomi in una variabile 
	 * @param C Vettore dei coefficienti
	 * 			(a partire dal meno significativo)
	 * @retun Vettore dei coefficienti di Bernstein
	 */
	public static Matrix toBernsteinMono(Matrix C) {
		
		// Grado del polinomio
		int degree = C.getColumnDimension()-1;
		
		// Matrice delle soluzioni parziali
		Matrix M = new Matrix(degree+1,degree+1);
		//double[][] M = new double[C.size()][C.size()];
		
		// Calcolo di termini in prima riga
		for(int j = 0; j <= degree; ++j) {
			M.set(0,j,C.get(0,j)/Binomial.binomial(degree,j));
		}
		
		int index = 1;
		
		// Calcolo i restanti termini con la formula ricorsiva
		for(int k = 1; k < degree + 1; ++k) {
			for(int j = 0; j < degree + 1 - index; ++j) {
				M.set(k,j,M.get(k-1, j+1)+M.get(k-1, j));
			}
			++index;
		}
		
		// Inizializzo il vettore risultato
		Matrix result = new Matrix(1,degree+1);
		
		/** Calcolo della soluzione a partire dalla 
		  matrice delle soluzioni parziali */
		for(int k = 0; k < degree + 1; ++k) {
			result.set(0,k,M.get(k,0));
		}
		
		return result;
	}
	
	
	/**********************************************************/
	/** Metodo per il ricalcolo dei coefficienti di Bernstein
	 * sull'intervallo dimezzato
	 * @param degree Grado del polinomio
	 * @param pBernstein Vettore dei coefficienti di Bernstein
	 * @return Vettore dei coefficienti ricalcolati sul nuovo intervallo
	 */
	public static Matrix NatarajRicalcolaBernstein(int degree, Matrix pBernstein) {
		
		// Definisco la matrice delle soluzioni parziali
		Matrix M = new Matrix(degree+1, degree+1);
		
		// Inserisco i primi coefficienti (noti) nella matrice
		for(int i = 0; i<= degree; ++i)
			M.set(0, i, pBernstein.get(0, i));
		
		for(int k = 1; k<=degree; ++k) {
			for(int i = 0; i<=degree; ++i) {
				// Applicazione Formula di Nataraj
				if(k > i)
					M.set(k, i, M.get(k-1, i));
				else // k <= i
					M.set(k, i, (0.5*M.get(k-1, i-1) + 0.5*M.get(k-1, i) ) );
			}
		}
		
		/** Matrice risultato: sulla prima riga i coefficienti ricalcolati sulla prima
		  metà dell'intervallo, mentre sulla seconda riga i coefficienti ricalcolati sulla
		  seconda metà dell'intervallo */
		Matrix result = new Matrix(2, degree+1);
		
		for(int i=0; i<=degree; ++i) {
			result.set(0, i, M.get(degree, i));
			result.set(1, i, M.get(degree-(i), degree));
		}
		
		return result;
	}
	
	
	/**********************************************************/
	/** Metodo per la restrizione di domini per vincoli in una
	 *  sola variabile.
	 * @param polinomio Matrice contenente i coefficienti in forma di potenze
	 * @param val Valore del vincolo
	 * @param toll Tolleranza asseganta
	 * @return risultato Matrice 2xN contenente gli intervalli accettati
	 */
	public static Matrix NatarajCheckConstraint(Matrix polinomio, double val, double toll) {
		// Porto in forma di Bernstein il polinomio
		Matrix pBernstein = toBernsteinMono(polinomio);
		
		// Calcolo il massimo e il minimo coefficiente di Bernstein
		double massimo = maximum(pBernstein);
		double minimo = minimum(pBernstein);
		
		// Definisco il vettore del dominio iniziale
		double[] dom = new double[2];
		dom[0] = 0.0;
		dom[1] = 1.0;
		
		if(massimo > val && minimo > val) {
			Matrix risultato = new Matrix(2,1);
			risultato.set(0, 0, 0.0);
			risultato.set(1, 0, 0.0);
			return risultato;
		} else {
			if(massimo < val && minimo < val) {
				Matrix risultato = new Matrix(2,1);
				risultato.set(0, 0, dom[0]);
				risultato.set(1, 0, dom[1]);
				return risultato;
			}
			else {
				Matrix risultato_p = 
						NatarajBiseziona(polinomio.getColumnDimension()-1,pBernstein, dom, val,toll);
				
				// Elimino gli zeri degli intervalli non accettati che ho inserito nella matrice
				int zeri_spuri = 0;
				for(int k = 0; k<risultato_p.getColumnDimension(); ++k)
					if(risultato_p.get(0,k) == 0 && risultato_p.get(1, k) == 0)
						zeri_spuri++;
				
				Matrix risultato = new Matrix(2,risultato_p.getColumnDimension()-zeri_spuri);
				int kk = 0;
				for(int k = 0; k<risultato_p.getColumnDimension(); ++k)
					if(risultato_p.get(0,k) != 0 || risultato_p.get(1, k) != 0) {
						risultato.set(0, kk, risultato_p.get(0, k));
						risultato.set(1, kk, risultato_p.get(1, k));
						++kk;
					}
				return risultato;
			}
		}
	}
	
	
	/**********************************************************/
	/** Function ausiliaria della NatarajCheckConstraint che compie la
	 * suddivisione del dominio
	 * @param degree Grado del polinomio
	 * @param pBernstein Vettore dei coefficienti di Bernstein
	 * @param dom Vettore del dominio
	 * @param val Valore che costituisce il vincolo
	 * @param toll Tolleranza assegnata
	 * @return Vettore di intervalli accettati
	 */
	private static Matrix 
		NatarajBiseziona(int degree, Matrix pBernstein, double[] dom, double val, double toll) {
		
		// Calcolo il massimo e il minimo coefficiente di Bernstein
		double massimo = maximum(pBernstein);
		double minimo = minimum(pBernstein);
		
		if(massimo > val && minimo > val) {
			Matrix risultato = new Matrix(2,1);
			risultato.set(0, 0, 0.0);
			risultato.set(1, 0, 0.0);
			return risultato;
		} else {
			if(massimo < val+toll && minimo < val+toll) {
				Matrix risultato = new Matrix(2,1);
				risultato.set(0, 0, dom[0]);
				risultato.set(1, 0, dom[1]);
				return risultato;
			}
			else {
				// Calcolo il punto medio
				double middle = (dom[0] + dom[1])/2;
				
				double[] dom1 = new double[2];
				double[] dom2 = new double[2];
				dom1[0] = dom[0];
				dom1[1] = middle;
				dom2[0] = middle;
				dom2[1] = dom[1];
				
				// Ricalcolo i coefficienti di Bernstein sui due sottointervalli
				Matrix ricalcolo = NatarajRicalcolaBernstein(degree, pBernstein);
				
				// Ricorro sulla prima metà
				Matrix res1 = 
					NatarajBiseziona(degree, ricalcolo.getMatrix(0, 0, 0, ricalcolo.getColumnDimension()-1), dom1, val, toll);
				
				// Ricorro sulla seconda metà
				Matrix res2 = 
					NatarajBiseziona(degree, ricalcolo.getMatrix(1, 1, 0, ricalcolo.getColumnDimension()-1), dom2, val, toll);
				
				// Unifico i risultati in una matrice unica
				Matrix final_res = new Matrix(2, res1.getColumnDimension()+res2.getColumnDimension());
				
				for(int i = 0; i<res1.getColumnDimension(); ++i) {
					final_res.set(0, i, res1.get(0, i));
					final_res.set(1, i, res1.get(1, i));
				}
				
				int k = 0;
				for(int i = res1.getColumnDimension();
						i<(res1.getColumnDimension()+res2.getColumnDimension()); ++i) {
					final_res.set(0, i, res2.get(0, k));
					final_res.set(1, i, res2.get(1, k));
					++k;
				}
				
				return final_res;
			}
		}
	}
	
	
	/************************************************************/
	/** Function che, data una Matrix restituisce l'elemento 
	 * massimo che essa contiene
	 * @param A Matrice da scansire
	 * @return Elemento massimo della matrice
	 */
	public static double maximum(Matrix A) {
		double res = A.get(0, 0);
		for(int i = 0; i <= A.getRowDimension()-1; ++i) {
			for(int j = 0; j <= A.getColumnDimension()-1; ++j) {
				if(res < A.get(i, j))
					res = A.get(i, j);
			}
		}
		return res;
	}
	
	
	/************************************************************/
	/** Function che, data una Matrix restituisce l'elemento 
	 * minimo che essa contiene
	 * @param A Matrice da scansire
	 * @return Elemento minimo della matrice 
	 */
	public static double minimum(Matrix A) {
		double res = A.get(0, 0);
		for(int i = 0; i <= A.getRowDimension()-1; ++i) {
			for(int j = 0; j <= A.getColumnDimension()-1; ++j) {
				if(res > A.get(i, j))
					res = A.get(i, j);
			}
		}
		return res;
	}
	
	
	/***********************************************************/
	/**  Metodo per la conversione da forma di potenze in forma 
	 * di Bernstein per polinomi in due variabili 
	 * @param C Matrice dei coefficienti del polinomio in forma di potenze
	 * @param dom Matrice dei domini delle due variabili
	 * @return Matrice dei coefficienti di Bernstein corrispondenti
	 */	
	public static Matrix toBernsteinDuo(Matrix C, Matrix dom) {
		
		// Vettore dei gradi assoluti delle due variabili
		int[] g_ass = new int[2];
		g_ass[0] = 0;
		g_ass[1] = 0;
		
		// Calcolo le diemensioni della matrice in input
		int[] dimM = new int[2];
		dimM[0] = C.getRowDimension();
		dimM[1] = C.getColumnDimension();
		
		// Calcolo i gradi assoluti delle singole variabili
		for(int i = 0; i<dimM[0]; ++i) {
			for(int j = 0; j<dimM[1]; ++j) {
				if(g_ass[0] < i && C.get(i, j) != 0)
		            g_ass[0] = i;
		        
		        if(g_ass[1] < j && C.get(i,j) != 0)
		            g_ass[1] = j;
		        
			}
		}
		
		/** Inizializzo la matrice finale dei coefficienti di Bernstein
		  del polinomio in due variabili */
		Matrix B_f = new Matrix(dimM[0], dimM[1]);
		
		for(int i = 0; i<dimM[0]; ++i) {
			for(int j = 0; j<dimM[1]; ++j) {
				/** Se l'elemento nella matrice è non nullo e non mi trovo nè sulla
		          prima riga né sulla prima colonna (non sono monomi insomma) */
				if(C.get(i,j) != 0 && (i != 0 || j != 0)) {
					Matrix B1 = monomio(i, g_ass[0], dom.get(0,0), dom.get(1,0));
					Matrix B2 = monomio(j, g_ass[1], dom.get(0,1), dom.get(1,1));
					//B1.print(2, 2);
					//B2.print(2, 2);
					Matrix B = (B1.transpose()).times((C.get(i,j))).times(B2);
					B_f = B_f.plus(B);
				}
				else {
					if(C.get(i,j) != 0 && i == 0 && j != 0) {
						Matrix B2 = monomio(j-1,g_ass[1],dom.get(0,1),dom.get(1,1));
						Matrix unit = new Matrix(g_ass[0],1,1);
			            B_f = B_f.plus( (unit.times(C.get(i,j))).times(B2) );
					}
					else {
						if(C.get(i,j) != 0 && j == 0 && i != 0) {
							Matrix B1 = monomio(i-1,g_ass[0],dom.get(0,0),dom.get(1,0));
							Matrix unit = new Matrix(1,g_ass[1],1);
				            B_f = B_f.plus( B1.times(C.get(i,j)).times(unit) );
						}
						else {
							if(C.get(i,j) != 0 && j == 0 && i == 0) {
								Matrix unit1 = new Matrix(g_ass[0]+1,1,1);
								Matrix unit2 = new Matrix(1,g_ass[1]+1,1);
								B_f = B_f.plus( unit1.times(C.get(1,1)).times(unit2) );
							}
						}
					}
				}
			} // fine primo for
		} // fine secondo for
		
		
		return B_f;
	}
	
	
	/*****************************************************************************/
	/** Metodo ausiliario della toBernsteinDuo per il calcolo dei
	 * coefficienti di Bernstein sui singoli monomi 
	 * @param g_rel Grado realativo 
	 * @param g_ass Grado assoluto
	 * @param minimo Estremo inferiore del dominio
	 * @param massimo Estremo superiore del dominio 
	 * @return Vettore dei coefficienti di Bernstein calcolato sul monomio uni-variabile
	 */
	public static Matrix monomio(int g_rel, int g_ass, double minimo, double massimo) {
		double[] b = new double[g_ass+1];
		
		if(g_rel == g_ass) { // PRIMO CASO
			for(int k = 0; k<=g_ass; ++k) {
				b[k] = Math.pow(minimo, g_rel - k)*Math.pow(massimo, k);
			}
		}
		else { // SECONDO CASO (g_ass > g_rel)
			int r = g_ass - g_rel;
			for(int k = 0; k <= g_ass; ++k) {
				int j = Math.max(0, k-r);
				int jj = Math.min(k, g_rel);
				b[k] = 0;
				
				for(int i = j; i<=jj; ++i) {
					double res1 = (Binomial.binomial(r, k-i)*Binomial.binomial(g_rel, i));
					double res2 = (Binomial.binomial(g_rel+r, k));
					double res3 = (Math.pow(minimo, g_rel - i)*Math.pow(massimo, i));
					b[k] = b[k] + (res1/res2)*res3;
				}
				
				
			}
		}

		Matrix result = new Matrix(b, 1);
		return result;
	}
	
	
	
	/***********************************************************************/
	/** Metodo per la riduzione di domini per vincoli
	 * a due dimensioni. 
	 * @param polinomio coefficienti del polinomio in forma di potenze
	 * @param dom matrice 2xN dei domini delle singole variabili
	 * @param val valore che costituisce il vincolo
	 * @param toll tolleranza assegnata
	 * @return matrice 2xN degli intervalli validi del dominio
	 * EURISTICA: taglio il dominio della variabile con dominio maggiore
	 */
	public static Matrix checkConstraintDuo(Matrix polinomio, Matrix dom1, Matrix dom2, double val, double toll) {
		
		/** Ricavo le dimensioni della matrice (che deve essere minima) */
		int[] dim = new int[2];
		dim[1] = polinomio.getRowDimension();
		dim[2] = polinomio.getColumnDimension();
		
		/** Porto la matrice di coefficienti in base di Bernstein */
		Matrix dom = new Matrix(2, 2);
		dom.set(0, 0, dom1.get(0, 0));
		dom.set(1, 0, dom1.get(1, 0));
		dom.set(0, 1, dom2.get(0, 0));
		dom.set(1, 1, dom2.get(1, 0));
		Matrix polB = toBernsteinDuo(polinomio, dom);
		
		/** Calcolo il massimo e il minimo coefficiente di Bernstein */
		double massimo = maximum(polB);
		double minimo = minimum(polB);
		
		if(massimo > val + toll && minimo > val + toll) {
		// Vincolo mai verificato
			Matrix risultato = new Matrix(4,1);
			risultato.set(0, 0, 0.0);
			risultato.set(1, 0, 0.0);
			risultato.set(2, 0, 0.0);
			risultato.set(3, 0, 0.0);
			return risultato;
		}
		else {
			if(massimo < val + toll && minimo < val + toll) {
			// Vincolo sempre verificato
				Matrix risultato = new Matrix(4,1);
				risultato.set(0, 0, dom.get(0, 0));
				risultato.set(1, 0, dom.get(1, 0));
				risultato.set(2, 0, dom.get(0, 1));
				risultato.set(3, 0, dom.get(1, 1));
				return risultato;
			}
			else { // Risultato non immediatamente deducibile
				Matrix D_res = NatarajBisezionaDuo(polB, dim[1] - 1, dim[2] - 1, dom1, dom2, val, toll);
			}
		}
		
		// TODO: da sostituire
		return polinomio;
	}
	
	
	
	/********************************************************/
	/** Metodo che permette di compiere la suddivisione
	 * del dominio per vincoli in due variabili.
	 * @param polB Matrix di coefficienti di Bernstein
	 * @param grado1 Grado massimo della prima variabile
	 * @param grado2 Grado massimo della seconda variabile
	 * @param dom1 Dominio della prima variabile
	 * @param dom2 Dominio della seconda variabile
	 * @param val Valore che costituisce il vincolo
	 * @param toll Tolleranza assegnata
	 */
	public static Matrix NatarajBisezionaDuo(Matrix polB, int grado1, int grado2, Matrix dom1, Matrix dom2, double val, double toll) {
		
		// Calcolo il massimo e minimo coefficiente di Bernstein
		double massimo = maximum(polB);
		double minimo = minimum(polB);
		
		if(massimo > val && minimo > val) {
			// Vincolo mai verificato
				Matrix risultato = new Matrix(4,1);
				risultato.set(0, 0, 0.0);
				risultato.set(1, 0, 0.0);
				risultato.set(2, 0, 0.0);
				risultato.set(3, 0, 0.0);
				return risultato;
			}
			else {
				if(massimo < val + toll && minimo < val + toll) {
				// Vincolo sempre verificato
					Matrix risultato = new Matrix(4,1);
					risultato.set(0, 0, dom1.get(0, 0));
					risultato.set(1, 0, dom1.get(1, 0));
					risultato.set(2, 0, dom2.get(0, 1));
					risultato.set(3, 0, dom2.get(1, 1));
					return risultato;
				}
				else {
					if(Math.abs(dom1.get(0, 1) - dom1.get(0, 0)) > (dom2.get(0, 1) - dom2.get(0, 0))) {
					// Taglio secondo la prima variabile
						double middle = (dom1.get(0,0) + dom1.get(1, 0))/2;
						Matrix domA = new Matrix(2,1);
						Matrix domB = new Matrix(2,1);
						domA.set(0, 0, dom1.get(0, 0));
						domA.set(0, 1, middle);
						domB.set(0, 0 ,middle);
						domB.set(0,1, dom1.get(1, 0));
						
						// Ricalcolo i coefficienti di Bernstein
						Matrix polBB = NatarajRicalcolaBernsteinDuo(polB, grado1, 1);
						Matrix polB1 = polBB.getMatrix(0, polB.getRowDimension()-1, 0, polB.getColumnDimension()-1);
						Matrix polB2 = polBB.getMatrix(0, polB.getRowDimension()-1, polB.getColumnDimension(), polBB.getColumnDimension()-1);
						
						// Continuo la suddivisione
						Matrix D1 = NatarajBisezionaDuo(polB1, grado1, grado2, domA, dom2, val, toll);
						Matrix D2 = NatarajBisezionaDuo(polB2, grado1, grado2, domB, dom2, val, toll);
						
						// Concateno i risultati ottenuti
						Matrix D_res = new Matrix(4,D1.getColumnDimension() + D2.getColumnDimension());
						
						for(int i = 0; i < D1.getColumnDimension(); ++i) {
							D_res.set(0, i, D1.get(0, i));
							D_res.set(1, i, D1.get(1, i));
							D_res.set(2, i, D1.get(2, i));
							D_res.set(3, i, D1.get(3, i));
						}
						
						int k = 0;
						for(int i = D1.getColumnDimension();
								i < (D1.getColumnDimension() + D2.getColumnDimension()); ++i) {
							D_res.set(0, i, D2.get(0, k));
							D_res.set(1, i, D2.get(1, k));
							D_res.set(2, i, D2.get(2, k));
							D_res.set(3, i, D2.get(3, k));
							++k;
						}
						
						return D_res;
					}
					else {
					// Taglio verso la seconda variabile
						double middle = (dom2.get(0,0) + dom2.get(1, 0))/2;
						Matrix domA = new Matrix(2,1);
						Matrix domB = new Matrix(2,1);
						domA.set(0, 0, dom2.get(0, 0));
						domA.set(0, 1, middle);
						domB.set(0, 0 ,middle);
						domB.set(0,1, dom2.get(1, 0));
						
						// Ricalcolo i coefficienti di Bernstein
						Matrix polBB = NatarajRicalcolaBernsteinDuo(polB, grado2, 2);
						Matrix polB1 = polBB.getMatrix(0, polB.getRowDimension()-1, 0, polB.getColumnDimension()-1);
						Matrix polB2 = polBB.getMatrix(0, polB.getRowDimension()-1, polB.getColumnDimension(), polBB.getColumnDimension()-1);
						
						// Continuo la suddivisione
						Matrix D1 = NatarajBisezionaDuo(polB1, grado1, grado2, dom1, domA, val, toll);
						Matrix D2 = NatarajBisezionaDuo(polB2, grado1, grado2, dom1, domB, val, toll);
						
						// Concateno i risultati ottenuti
						Matrix D_res = new Matrix(4, D1.getColumnDimension() + D2.getColumnDimension());
						
						for(int i = 0; i < D1.getColumnDimension(); ++i) {
							D_res.set(0, i, D1.get(0, i));
							D_res.set(1, i, D1.get(1, i));
							D_res.set(2, i, D1.get(2, i));
							D_res.set(3, i, D1.get(3, i));
						}
						
						int k = 0;
						for(int i = D1.getColumnDimension();
								i < (D1.getColumnDimension() + D2.getColumnDimension()); ++i) {
							D_res.set(0, i, D2.get(0, k));
							D_res.set(1, i, D2.get(1, k));
							D_res.set(2, i, D2.get(2, k));
							D_res.set(3, i, D2.get(3, k));
							++k;
						}
						
						return D_res;
					}
				}
			}
	}
	
	
	/***********************************************************/
	/** Metodo per il ricalcolo dei coefficienti di Bernstein
	 * di polinomi in due variabili
	 * @param polB Matrice dei coefficienti di Bernstein
	 * @param degree Grado massimo della variabile di taglio
	 * @param taglio Direzione del taglio
	 * @return
	 */
	public static Matrix NatarajRicalcolaBernsteinDuo(Matrix polB, int degree, int taglio) {
		// Occorre una matrice tridimensionale per applicare agilmente la formula
		
		// TODO: da modificare
		return polB;
	}
	
	
	
}
