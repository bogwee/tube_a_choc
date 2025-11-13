# Simulation numérique du tube à choc (Projet MAM4 - Sujet 8)

<br>

## Objectif du projet

Ce projet a pour but de simuler le comportement d’un gaz dans un **tube à choc** après la rupture d’un diaphragme séparant deux zones de pression différentes.  
Il s’agit d’un **problème de Riemann** typique, dont la résolution fait intervenir les **équations d’Euler** et la gestion de **discontinuités physiques** : onde de choc, onde de détente et ligne de contact.

L’objectif principal est de développer une **simulation numérique en C++** pour le cas **1D**, puis d’étendre la modélisation au **cas 2D**.  
La partie visualisation des résultats est réalisée en **Python**, afin d’afficher les profils de densité, pression et vitesse au cours du temps.

<br>

## Contexte physique

Le tube à choc est un dispositif utilisé pour étudier les écoulements compressibles.  
Il est constitué d’un tube fermé, divisé en deux parties par un diaphragme :

- **Partie gauche (gaz haute pression)** : densité initiale $\rho_L$, pression $p_L$, température $T_L$  
- **Partie droite (gaz basse pression)** : densité $\rho_R$, pression $p_R < p_L$, température $T_R$

Avant la rupture, le fluide est au repos ($U = 0$).  
Au moment où le diaphragme cède ($t = 0$), trois phénomènes se produisent :

1. Une **onde de détente** se propage vers la gauche (dans le gaz à haute pression).  
2. Une **onde de choc** se propage vers la droite (dans le gaz à basse pression).  
3. Une **discontinuité de contact** sépare les deux zones intermédiaires.

Ces phénomènes constituent un exemple fondamental d’écoulement compressible non linéaire et sont étudiés pour valider les méthodes numériques utilisées en mécanique des fluides.

<br>

## Schéma du tube à choc

*(Espace réservé pour insérer ton image du tube à choc :)*

![Tube à choc – insérer ici ton image](images/tube.png)

<br>

## Modélisation mathématique

Le comportement du gaz est régi par les **équations d’Euler** pour un fluide compressible sans viscosité :

$$
\frac{\partial W}{\partial t} + \frac{\partial F(W)}{\partial x} = 0
$$

où :

$$
W =
\begin{pmatrix}
\rho \\
\rho U \\
E
\end{pmatrix},
\qquad
F(W) =
\begin{pmatrix}
\rho U \\
\rho U^2 + p \\
(E + p)U
\end{pmatrix}
$$

Les variables représentent :
- $\rho$ : densité du gaz  
- $U$ : vitesse moyenne  
- $E$ : énergie totale par unité de volume  
- $p$ : pression  

La pression suit la **loi des gaz parfaits** :

$$
p = \rho R T
$$

et l’énergie totale s’écrit :

$$
E = \frac{p}{\gamma - 1} + \frac{1}{2} \rho U^2
$$

où $\gamma$ est le rapport des capacités calorifiques ($\gamma = 1.4$ pour l’air).

<br>

## Conditions initiales (cas 1D)

Le tube est divisé en deux zones initialement distinctes par un diaphragme placé en $x_0 = 0.5$ :

$$
W(x,0) =
\begin{cases}
(\rho_L, \rho_L U_L, E_L), & x < x_0 \\
(\rho_R, \rho_R U_R, E_R), & x > x_0
\end{cases}
$$

Dans le **cas de Sod**, les valeurs typiques sont :

| Zone | Densité $\rho$ | Pression $p$ | Vitesse $U$ | Position |
|------|----------------|---------------|--------------|-----------|
| Gauche (L) | 8 | $10/\gamma$ | 0 | $x < 0.5$ |
| Droite (R) | 1 | $1/\gamma$ | 0 | $x > 0.5$ |

Le domaine spatial est $[0,1]$ et on observe la solution à $t = 0.2$.

<br>

## Discrétisation numérique

Le domaine spatial est discrétisé en $N_x$ mailles, de taille :

$$
\Delta x = \frac{1}{N_x - 1}
$$

Le pas de temps est choisi selon le **critère de stabilité CFL** :

$$
\Delta t = \text{CFL} \cdot \frac{\Delta x}{|U| + a}
$$

avec $a = \sqrt{\gamma p / \rho}$ la célérité du son.

Ce critère permet d’assurer que les informations physiques ne se propagent pas plus vite que le maillage numérique ne peut les suivre.

<br>

## Méthodes de résolution

Plusieurs schémas numériques peuvent être employés pour résoudre les équations d’Euler.



### 1. Schémas centrés

Ces méthodes utilisent des approximations symétriques pour estimer les flux.  
Elles offrent une bonne précision pour les régions où les variables sont continues, mais peuvent produire des oscillations près des discontinuités.  

Exemples :
- **Lax-Wendroff**
- **MacCormack**

Une **dissipation artificielle** peut être ajoutée pour stabiliser la solution.



### 2. Schéma décentré (upwind)

Ce type de schéma tient compte de la direction de propagation des ondes.  
Il est plus stable et mieux adapté aux chocs.  

Un exemple typique est le **schéma de Roe**, basé sur une **approximation de Riemann linéarisée**.

Le flux numérique de Roe s’écrit :

$$
F(W_L, W_R) = \frac{1}{2}\left[F(W_L) + F(W_R) - |A_{L,R}| (W_R - W_L)\right]
$$

où la matrice de Roe $A_{L,R}$ dépend des valeurs moyennes :

$$
U_{L,R} = \frac{\sqrt{\rho_L}U_L + \sqrt{\rho_R}U_R}{\sqrt{\rho_L} + \sqrt{\rho_R}}
$$

$$
H_{L,R} = \frac{\sqrt{\rho_L}H_L + \sqrt{\rho_R}H_R}{\sqrt{\rho_L} + \sqrt{\rho_R}}
$$

avec $H = (E + p) / \rho$ l’enthalpie totale.

<br>

## Résultats attendus

La simulation doit permettre d’observer :

- La formation d’une onde de détente à gauche,  
- La propagation d’un choc à droite,  
- La zone intermédiaire séparée par une discontinuité de contact.

Les variables d’intérêt sont la densité $\rho(x, t)$, la pression $p(x, t)$ et la vitesse $U(x, t)$.

Les résultats seront sauvegardés dans un **fichier CSV** sous la forme de **vecteurs** représentant les valeurs discrètes de chaque variable physique.  
Ces fichiers seront ensuite importés et manipulés en **Python** pour la visualisation et l’analyse.

<br>

## Visualisation des résultats

La visualisation est réalisée à l’aide de **Python**, en utilisant les bibliothèques `matplotlib`, `numpy` et `pandas`.

Les scripts permettent de :
- Lire les données exportées depuis le programme C++ (vecteurs enregistrés dans un fichier `.csv`),
- Tracer les profils de densité $\rho(x)$, vitesse $U(x)$ et pression $p(x)$,
- Comparer les résultats numériques avec la solution analytique.

Exemple de visualisation (à insérer) :

![Résultats de la simulation](images/resultats.png)

<br>

## Extension vers le cas 2D

L’étape suivante consiste à étendre le modèle pour intégrer la dimension temporelle et spatiale complète, ce qui permettra de simuler des phénomènes tels que :

- La propagation d’ondes de choc en deux dimensions,  
- Les interactions entre plusieurs chocs,  
- L’étude de la stabilité et de la précision des schémas numériques en 2D.

Cette extension nécessitera une généralisation du flux bidimensionnel et une gestion plus fine des conditions aux limites.

<br>

## Bibliothèques utilisées

### En C++ :
- `cmath` : opérations mathématiques et fonctions standards,  
- `vector` : stockage dynamique des données,  
- `fstream` : lecture et écriture des fichiers CSV.

### En Python :
- `matplotlib` : visualisation des résultats,  
- `numpy` : manipulation de tableaux numériques,  
- `pandas` : lecture et traitement des fichiers CSV.

<br>

## Équipe projet

Projet réalisé dans le cadre du module **MAM4 - Calcul scientifique**  
Université Côte d’Azur  

Sujet 8 : *Étude du tube à choc – Problème de Riemann et discontinuités*  

Langages utilisés :
- **C++** : simulation et calcul numérique  
- **Python** : visualisation et analyse graphique
