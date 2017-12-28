#!/bin/bash
#da correggere!


/home/nicholas/tesi/code/get-film-screenshots action /home/nicholas/tesi/input/train/screenshots
/home/nicholas/tesi/code/get-film-screenshots drama /home/nicholas/tesi/input/train/screenshots
/home/nicholas/tesi/code/get-film-screenshots list-comedy /home/nicholas/tesi/input/train/screenshots
/home/nicholas/tesi/code/get-film-screenshots animation /home/nicholas/tesi/input/train/screenshots
/home/nicholas/tesi/code/get-film-screenshots action /home/nicholas/tesi/input/test/screenshots
/home/nicholas/tesi/code/get-film-screenshots drama /home/nicholas/tesi/input/test/screenshots
/home/nicholas/tesi/code/get-film-screenshots list-comedy /home/nicholas/tesi/input/test/screenshots
/home/nicholas/tesi/code/get-film-screenshots animation /home/nicholas/tesi/input/test/screenshots
python3 /home/nicholas/tesi/code/python/create_lmdb.py
compute_image_mean -backend=lmdb /home/nicholas/tesi/input/train_lmdb /home/nicholas/tesi/input/mean.binaryproto
caffe train -solver /home/nicholas/tesi/caffe-model/solver_1.prototext 2>&1 | tee /home/nicholas/tesi/caffe-model/model_train.log
