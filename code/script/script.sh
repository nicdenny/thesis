#!/bin/bash

#extracting screenshot from trailers for train data set
/thesis/code/script/get-film-screenshots action /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots drama /thesis/input/train/screenshots
/thesis/code/script/get-film-screenshots list-comedy /thesis/input/train/screenshots
/code/script/get-film-screenshots animation /thesis/input/train/screenshots

#extracting screenshot from trailers for validation data set
/thesis/code/script/get-film-screenshots action /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots drama /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots list-comedy /thesis/input/test/screenshots
/thesis/code/script/get-film-screenshots animation /thesis/input/test/screenshots


#python3 /home/nicholas/tesi/code/python/create_lmdb.py
#compute_image_mean -backend=lmdb /home/nicholas/tesi/input/train_lmdb /home/nicholas/tesi/input/mean.binaryproto
#caffe train -solver /home/nicholas/tesi/caffe-model/solver_1.prototext 2>&1 | tee /home/nicholas/tesi/caffe-model/model_train.log
