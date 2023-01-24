from sklearn import datasets
from sklearn.model_selection import train_test_split
import numpy as np


def sigmoid(val):
    return 1 / (1 + np.exp(-val))


def SigmoidDerivative(x):
    return np.multiply(x, 1-x)


class Network:
    # Create Directed Acyclic Network of given number layers.
    def __init__(self, num_layers):
        np.random.seed(10)
        # Hidden layer - number of columns is the number of hidden layers - the input layer
        self.hidden = np.random.normal(scale=0.5, size=(4, num_layers))
        # Output layer
        self.out = np.random.normal(scale=0.5, size=(num_layers, 3))

    def forward(self, train):
        # initialize each weight with the values min_value=-0.5, max_value=0.5,
        # The activation function used is sigmoid function
        self.hidden_layer = sigmoid(np.dot(train, self.hidden))
        self.out_layer = sigmoid(np.dot(self.hidden_layer, self.out))
        return self.out_layer


def BackPropagationLearner(X, y, net, learning_rate, epochs, error_print=False):
    # labels has shape (150,) while predictions has shape (150,3)
    # label array needs to be arranged
    arranged_y = np.zeros((y.size, 3))             # there are 3 labels
    arranged_y[np.arange(y.size), y] = 1

    for epoch in range(epochs):
        # Forward pass
        net.forward(X)

        # Error for the MSE cost function
        if error_print:
            error = ((net.out_layer - arranged_y) ** 2).sum() / (2 * net.out_layer.size)
            print("Mean squared error is: {:.4f}".format(error))

        # Backward pass
        d_out = np.dot(net.hidden_layer.T, (arranged_y - net.out_layer) * SigmoidDerivative(net.out_layer))
        d_hidden = np.dot(X.T, np.dot((arranged_y - net.out_layer) * SigmoidDerivative(net.out_layer),
                          net.out.T) * SigmoidDerivative(net.hidden_layer))
        #  Update weights
        net.hidden += d_hidden * learning_rate
        net.out += d_out * learning_rate

    error = ((net.out_layer - arranged_y) ** 2).sum() / (2 * net.out_layer.size)
    print("Mean squared error is: {:.4f}".format(error))
    return net


def NeuralNetLearner(X, y, hidden_layer_sizes=None, learning_rate=0.01, epochs=100):
    # construct a raw network and call BackPropagationLearner
    if hidden_layer_sizes is None:
        hidden_layer_sizes = 3
    else:
        hidden_layer_sizes = len(hidden_layer_sizes)
    network = Network(hidden_layer_sizes)
    network = BackPropagationLearner(X, y, network, learning_rate, epochs)

    def predict(example):
        # forward pass
        prediction = network.forward(example)
        # find the max node from output nodes
        prediction = prediction.argmax(axis=0)
        return prediction

    return predict


if __name__ == '__main__':
    iris = datasets.load_iris()
    X = iris.data
    y = iris.target

    nNL = NeuralNetLearner(X, y)
    print(nNL([4.6, 3.1, 1.5, 0.2]))  # 0
    print(nNL([6.5, 3., 5.2, 2.]))  # 2


    """     #### To find the best configuration ####
    for hidden_layer_size in (2, 3, 4):
        for learning_rate in (0.1, 0.01, 0.001):
            print("Configuration layer size: {} - lr: {}".format(hidden_layer_size, learning_rate))
            NeuralNetLearner(X, y, hidden_layer_size, learning_rate)
    ## best configuration found as -> layer_size: 4, learning_rate: 0.1
    """