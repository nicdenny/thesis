#!/bin/bash

#extracting screenshot from trailers for train data set
/thesis/code/script/get-film-screenshots action /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots adventure /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots animation /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots comedy /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots crime /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots drama /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots fantasy /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots horror /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots romance /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots superhero /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots thriller /thesis/input/train/screenshots

#extracting screenshot from trailers for validation data set
/thesis/code/script/get-film-screenshots action /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots adventure /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots animation /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots comedy /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots crime /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots drama /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots fantasy /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots horror /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots romance /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots superhero /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots thriller /thesis/input/test/screenshots

/thesis/code/bezier/bezier_extractor /thesis/input/train/screenshots/*.png
/thesis/code/bezier/bezier_extractor /thesis/input/test/screenshots/*.png

echo "Moving bezier images...(train)"
mv /thesis/input/train/screenshots/*.jpg /thesis/input/train/bezier-image
echo "Moving bezier data...(train)"
mv /thesis/input/train/screenshots/*.txt /thesis/input/train/bezier-data

echo "Moving bezier images...(test)"
mv /thesis/input/train/screenshots/*.jpg /thesis/input/test/bezier-image
echo "Moving bezier data...(test)"
mv /thesis/input/train/screenshots/*.txt /thesis/input/test/bezier-data

#python3 /home/nicholas/tesi/code/python/create_lmdb.py
#compute_image_mean -backend=lmdb /home/nicholas/tesi/input/train_lmdb /home/nicholas/tesi/input/mean.binaryproto
#caffe train -solver /home/nicholas/tesi/caffe-model/solver_1.prototext 2>&1 | tee /home/nicholas/tesi/caffe-model/model_train.log
