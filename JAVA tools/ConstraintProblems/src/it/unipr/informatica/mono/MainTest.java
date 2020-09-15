package it.unipr.informatica.mono;

import java.util.ArrayList;
import java.util.Vector;

import Jama.*;

public class MainTest {

	public static void main(String[] args) {
		
		/********* Test sulla toBernsteinMono **********/
		
		Matrix prova = new Matrix(1,4);
		prova.set(0, 0, -8);
		prova.set(0, 1, 65);
		prova.set(0, 2, -150);
		prova.set(0, 3, 90);
		System.out.println("Dato il polinomio di coefficienti (a partire dal meno significativo):");
		prova.print(2, 0);
		
		Matrix result = BernsteinTools.toBernsteinMono(prova);
		System.out.println("Si ottiene la sua forma di Bernstein:");
		result.print(2, 3);
		
		/********* Test sulla NatarajRicalcolaBernstein **********/
		
		Matrix ricalcolo = BernsteinTools.NatarajRicalcolaBernstein(3, result);
		System.out.println("Spezzando l'intervallo [0,1], su cui il polinomio è definito, "
				+ "nei due intervalli [0,0.5] e [0.5,1]" + "\n"
				+ "si ottengono i seguenti nuovi coefficienti di Bernstein:");
		Matrix primaMeta = ricalcolo.getMatrix(0, 0, 0, ricalcolo.getColumnDimension()-1);
		Matrix secondaMeta = ricalcolo.getMatrix(1, 1, 0, ricalcolo.getColumnDimension()-1);
		primaMeta.print(2, 4);
		System.out.println("per il primo intervallo, mentre:");
		secondaMeta.print(2, 4);
		System.out.println("per il secondo intervallo.");
		
		/********* Test sulla NatarajCheckConstraint **********/
		
		System.out.println("Proviamo ora a ridurre il dominio del vincolo:\n");
		System.out.println("-8 + 65x -150x^2 + 90x^3 < 0 " + "per x appartenente a [0,1].\n");
		System.out.println("Si ottengono gli intervalli: ");
		Matrix vincolo = BernsteinTools.NatarajCheckConstraint(prova, 0, 0.3);
		vincolo.print(2, 4);
		
		/********** Test sulla toBernsteinDuo *****************/
		
		Matrix dueVar = new Matrix(4,3);
		dueVar.set(0, 0, 0);
		dueVar.set(0, 1, 0);
		dueVar.set(0, 2, 0);
		dueVar.set(1, 0, 0);
		dueVar.set(1, 1, -6);
		dueVar.set(1, 2, 0);
		dueVar.set(2, 0, 0);
		dueVar.set(2, 1, 0);
		dueVar.set(2, 2, 0);
		dueVar.set(3, 0, 0);
		dueVar.set(3, 1, 0);
		dueVar.set(3, 2, 1);
		
		Matrix dominion = new Matrix(2, 2);
		dominion.set(0, 0, 0);
		dominion.set(0, 1, 0);
		dominion.set(1, 0, 1);
		dominion.set(1, 1, 1);
		
		Matrix dueVarBernstein = BernsteinTools.toBernsteinDuo(dueVar, dominion);
		
		System.out.println("Testiamo ora la toBernsteinDuo. Consideriamo la matrice " 
				+ "rappresentante i coefficienti del polinomio in due variabili:");
		dueVar.print(2, 0);
		System.out.println("Con le variabili aventi rispettivamente dominio:");
		dominion.print(2, 0);
		System.out.println("Allora la sua forma di Bernstein è:");
		dueVarBernstein.print(2, 4);
		
		
		return;
	}

}
