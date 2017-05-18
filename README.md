# Reconstruction d'images par approche TV-TV^2
Les programmes implémentés sont décrits dans la thèse de K. Papafitsoros disponible
[ici](https://www.repository.cam.ac.uk/handle/1810/246692).

## Utilisation des programmes
### Débruitage d'images
Pour compiler ce programme, on utilise la commande

    make denoiser

Pour exécuter ce programme, on utilise

    bin/denoiser image l1 l2 a b

Le programme recherche le fichier image dans `img/input/`
Pour les paramètres par défaut, on peut utiliser

    bin/denoiser < input_denoiser

### Inpainting
POur compiler ce programme, on utilise la commande

    make inpainter

Pour exécuter ce programme, on utilise

    bin/inpainter image l1 l2 a b

Le programme recherche le fichier image dans `img/input/`

Pour les paramètres par défaut, on peut utiliser

    bin/inpainter < input_inpainter


### En cas d'erreur à la compilation
Créer les dossiers `ob/` et `bin/`.


## Futures versions
* Réorganisation du code
* Ordre de remplissage des piexels pour l'inpainting
