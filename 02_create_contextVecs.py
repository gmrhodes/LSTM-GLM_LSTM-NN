# This file creates the window-specific context vectors for the simulations.

#################################################### USER SET-UP ####################################################
#Survival model ("aft" or "cox")
surv_model = "aft"

#Prediction time 
t = 5

#Number of iterations
iters = 500

#Number of cores to parallelize across
num_cpu = 25

#Directory path of folder containing data generated by 01_generate_data.R
##Must contain sub-directories "aft" and "cox"
input_path = "D:/Simulations/data/"

#Directory path of folder to write KerasTuner data to
##Must contain sub-directories "aft" and "cox"
tuner_path = "D:/Simulations/tuner/"

#Directory path of output folder to write results to
##Must contain sub-directories "aft" and "cox"
output_path = "D:/Simulations/context_vectors/"


#################################################### SET-UP ####################################################
seed_number=5678
import os
os.environ['PYTHONHASHSEED']='0'
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
import random
random.seed(seed_number)
import numpy as np
np.random.seed(seed_number)
import kerastuner as kt
import tensorflow as tf
from tensorflow.keras import layers as tfkl
tf.random.set_seed(seed_number)
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
import pandas as pd
from sklearn import preprocessing
import multiprocessing as mp

#Padding value for missing time steps
padVal = -9999 

#LSTM autoencoder hyperparameters
epochs = 500 #training epochs
cv_dimension = 5 #context vector dimension


#################################################### FUNCTION TO SCALE COVARIATES ####################################################
def minMaxScale(df):
    #Select continuous longitudinal covariate
    contCovDF = df[["long_cov"]].copy()

    #Min-max scale data
    scaler = preprocessing.MinMaxScaler()
    dfStand = scaler.fit_transform(contCovDF)
    dfStand = pd.DataFrame(dfStand, columns=contCovDF.columns)
    
    #Replace scaled variable in original dataset
    df[["long_cov"]] = dfStand.values
    return df


#################################################### FUNCTION TO CREATE 3D COVARIATE ARRAY ####################################################
def covArray_construct(covDF):    
    #Create padded 3D numpy arrays of covariates for Keras
    gb = covDF.groupby(['id'])
    maxBlocks = gb['id'].size().max() 
    covArr = np.array([np.pad(frame['block'].values, pad_width=(0, maxBlocks-len(frame)), mode='constant', constant_values=padVal) for _,frame in gb]).reshape(-1, maxBlocks, 1)
    for col in covDF.columns[2:]:
        newArr = np.array([np.pad(frame[col].values, pad_width=(0, maxBlocks-len(frame)), mode='constant', constant_values=padVal) for _,frame in gb]).reshape(-1, maxBlocks, 1)    
        covArr = np.dstack((covArr,newArr))
    return covArr


#################################################### FUNCTIONS TO DEFINE LSTM AUTOENCODER #################################################### 
#Class returns the output of an LSTM layer stacked on a RepeatVector layer with the mask propogated through both layers
class lstm_bottleneck(tf.keras.layers.Layer):
    def __init__(self, lstm_units, time_steps, **kwargs):
        self.lstm_units = lstm_units
        self.time_steps = time_steps
        self.lstm_layer = tfkl.LSTM(units=lstm_units, activation='relu', return_sequences=False)
        self.repeat_layer = tfkl.RepeatVector(time_steps)
        super(lstm_bottleneck, self).__init__(**kwargs)

    def call(self, inputs):
        return self.repeat_layer(self.lstm_layer(inputs))

    def compute_mask(self, inputs, mask=None):
        return mask

#Class creates HyperModel object for LSTM autoencoder
class lstmHyperModel(kt.HyperModel):
    def __init__(self, padVal, timeSteps):
        self.padVal = padVal
        self.timeSteps = timeSteps
        
    def build(self, hp):
        input_layer = tfkl.Input(shape=(self.timeSteps, 1))
        x = tfkl.Masking(mask_value=self.padVal)(input_layer)
        x = lstm_bottleneck(lstm_units=cv_dimension, time_steps=self.timeSteps)(x)
        x = tfkl.LSTM(units=cv_dimension, activation='relu', return_sequences=True)(x)
        x = tfkl.TimeDistributed(tfkl.Dense(1))(x)
        lstm_ae = tf.keras.models.Model(inputs=input_layer, outputs=x)
        lstm_ae.compile(optimizer=tf.keras.optimizers.Adam(hp.Choice('learning_rate', values=[1e-3,1e-4])), loss='mse')
        return lstm_ae
    

#################################################### FUNCTION TO CONSTRUCT WINDOW-SPECIFIC CONTEXT VECTORS #################################################### 
#Function to construct window-specific context vectors for simulation ind
def context_vector_construct(longCov_lstm, longName, ind):
    #Set-up tensorflow session
    session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
    sess = tf.compat.v1.Session(config=session_conf)
    tf.compat.v1.keras.backend.set_session(sess)

    #Construct LSTM HyperModel & conduct cross-validation to select learning rate
    lstm_hyperMod = lstmHyperModel(padVal=padVal, timeSteps=longCov_lstm.shape[1])
    tuner = kt.tuners.RandomSearch(lstm_hyperMod, objective='val_loss', max_trials=2, seed=8, overwrite=True, 
                                           directory='{}{}/'.format(tuner_path, surv_model), project_name='iter{}'.format(ind))
    tuner.search(longCov_lstm, longCov_lstm, epochs=epochs, validation_split=0.2, verbose=0)
    
    #Retrieve model with "best" learning rate & fit to data
    lstm_best = tuner.get_best_models(num_models=1)[0]
    lstm_best.fit(longCov_lstm, longCov_lstm, epochs=epochs, verbose=0)

    #Extract context vectors from LSTM autoencoder 
    contextVec_model = tf.keras.Model(inputs=lstm_best.inputs, outputs=lstm_best.layers[2].output)
    cv = contextVec_model.predict(longCov_lstm)[:,0,:]  
    
    #Save variable names 
    contVecNames = list()
    for j in np.arange(1, cv_dimension+1, 1):
        contVecNames.append('{}_{}'.format(longName,j))
    
    #Reset tensorflow session & seeds 
    del tuner
    del lstm_hyperMod
    del lstm_best
    del contextVec_model
    tf.compat.v1.keras.backend.clear_session()
    tf.compat.v1.reset_default_graph()    
    random.seed(seed_number)
    np.random.seed(seed_number)
    tf.random.set_seed(seed_number) 
    
    #Return window-specific context vectors & associated variable names
    return (cv, contVecNames)      
      
    
################################################## FUNCTION TO PERFORM A SINGLE ITERATION #################################################### 
#Function to output window-specific context vectors for simulation ind
def iterFN(ind, timeDF, covDF):    
    random.seed(seed_number)
    np.random.seed(seed_number)
    tf.random.set_seed(seed_number)
    
    #Merge survival time and covariate data 
    simDF = timeDF.merge(covDF, on='id')
    
    #Min-max scale longitudinal covariate
    simDF = minMaxScale(simDF)
    
    #Keep only measurements taken prior to t on patients at risk at t
    simDF_t = simDF[(simDF.meas_time<t) & (simDF.Y>t)].copy()
    
    #Create response dataframe 
    respDF = simDF_t[["id", "base_cov", "delta_star", "Y"]].copy()
    respDF['restricted_residual_life'] = respDF['Y']-t
    respDF = respDF.drop_duplicates(subset='id', keep='first')
    respDF.index = range(0,len(respDF))
    
    #Create 3D covariate array (necessary for Keras)
    covDF = simDF_t[["id", "block", "long_cov"]].copy()
    covArr = covArray_construct(covDF)
    
    #Create context vector for longitudinal covariate
    longArr = np.delete(covArr,[0],2)
    longName = list(covDF.columns)[2]
    result = context_vector_construct(longArr, longName, ind)  
                
    #Combine baseline covariates, context vector, and response variables
    regressMatrix = pd.concat([respDF, pd.DataFrame(result[0], columns=result[1])], axis=1)

    #Return index & dataframe
    return(ind, regressMatrix)
    
    
#################################################### FUNCTION TO COLLECT PARALLELIZED RESULTS ####################################################
def get_iter_results(result):
    #Save index
    global index_list
    index_list.append(result[0])
    #Save data matrix
    global regressMat_list
    regressMat_list.append(result[1])
    
    
#################################################### MAIN ####################################################
if __name__ == '__main__':    
    #Read covariate data
    covDF = pd.read_csv('{}covariates.csv'.format(input_path)) 
    
    #Read survial time data for each simulation
    timeDF_list = list()
    for i in np.arange(1,iters+1):        
        timeDF = pd.read_csv('{}{}/survTimes{}.csv'.format(input_path,surv_model,i)) 
        timeDF_list.append(timeDF)        
    
    #Create lists to save context vector data sets & corresponding simulation indices
    index_list = list()
    regressMat_list = list()    
    
    #For each iteration, construct & output the window-specific context vectors
    pool = mp.Pool(num_cpu)
    asyncResults = [pool.apply_async(iterFN, args=(i+1, timeDF_list[i], covDF), callback=get_iter_results) for i in np.arange(0,iters)]       
    pool.close()
    pool.join()   
    
    #Write context vector data set for each simulation
    for i in np.arange(0,iters):        
        regressMat_list[i].to_csv('{}{}/contextVec{}.csv'.format(output_path,surv_model,index_list[i]), index=False)
    
    
        
        
    
    



   









