![Project Logo](assets/banner1.png)

![Coverage Status](assets/github.png)
![Coverage Status](assets/python.png)

<h1 align="center">
StereoChem : Explore the world of isomers through a fun game 
</h1>

<br>

### Collaborators: Catherina Carrer, Lucie Frenot, Mridhula Jayasankar ＆ Rania Doukkali
#### Practical Programming in Chemistry @ EPFL

## 📖 Package description 

Welcome to the world of isomers! If you are struggling to understand the chapter of isomers in chemistry, then this educational game is for you. 

Through an exciting interactive game, this package unables the user to get familiar with optical stereoisomers and additionally the chirality of a given molecule by making the learning enjoyable and engaging. Whether you are trying to understand this complicated chapter in chemistry or just wanting to explore a new aspect of this subject, this package offers a lively way to explore enantiomers, non-superposable mirror images, and diastereoisomers, non-mirror image optical isomers. 


## ⚛ What are isomers and chirality ? 

All chemical molecules can be represented in several ways one of which is the chemical formula, a universal method of notation. When two or more molecules have the same chemical formula but differ in the arrangement of their atoms giving rise to different physical and chemical properties, they are known as isomers.

Being linked to enantiomers, chirality represents the property of a molecule making it non superposable to its mirror image. Having a hard time finding them? It is easy find a carbon atom with four different constituents and you are all done. 


## 👩‍💻 Installation

Have we motivated you to come along the wonderful experience of StereoChem with us ? Start by creating a new environment! You may also give the environment a different name. 

```
conda create -n stereochem python=3.10 
```

```
conda activate stereochem
(conda_env) $ pip install .
```

If you need jupyter lab, install it 

```
(stereochem) $ pip install jupyterlab
```

As this package requires a few external packages please install the following if they don’t appear in your environment. In order to check if *rdkit3*, *streamlit*, *streamlit_ketcher* and *pubchempy* are already 
installed run:

```
(stereochem) $ conda list
```

If they don't appear, please install them individually: 

```
(stereochem) $ pip install rdkit 
(stereochem) $ pip install streamlit 
(stereochem) $ pip install streamlit_ketcher
(stereochem) $ pip install pubchempy 
```

## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:mridhula28/StereoChem`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:mridhula28/StereoChem.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(stereochem) $ pip install -e ".[test,doc]"
```

## 🔥 Usage 

Once all the boring installation work is done, the fun can start!
Open the `src` folder and open the `stereochem_project` file. Go to the terminal at type in: 

```
streamlit run StereoChem/src/stereochem/stereochem_project.py     
```

Executing this code will automatically bring you to the our web page. 
Here is a short user manual to help you get familiar with the web interface.

Start by entering the name of a molecule either by typing it in the sidebar and pressing on `Submit Name` or by drawing it directly on the interface and then clicking on `Apply` and `Submit drawing`. An example with ** butan-2-ol **

<img src="StereoChem/assets/name.png" alt="Alt text" width="500" />

Switch tabs and go to the ** Draw isomers ** tab. 
Guess all the stereoisomers on the interface and press on `Apply` and `Submit drawing` to validate your guess. If you are blocked press on `hint` or if you wish to give up and have the answers press on `Show Answers`

For each correct stereoisomer guessed, a small pop up card will appear. Try and the guess the iupac name of the molecule you just drew! 

<img src="StereoChem/assets/cards.png" alt="Alt text" width="500" />

Switch tabs and go to the ** Chirality ** tab. 
Click on the numbered boxes of the atoms you think are chiral.

<img src="StereoChem/assets/chirality.png" alt="Alt text" width="500" />

rien ne marche AHHHHHHHH

### 📚 Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox

