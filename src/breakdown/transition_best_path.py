import networkx as nx
import heapq

class TransitionBestPath:
    def __init__(self, stateTransitions, startState, maxDepthBest, maxDepthPossible):
        """Initialize TransitionBestPath object.

        Args:
            stateTransitions (dict): Dictionary representing state transitions.
            startState (str): Starting state of the Transition process.
        """
        self.stateTransitions = stateTransitions
        self.startState = startState
        self.maxDepthBest = maxDepthBest
        self.maxDepthPossible = maxDepthPossible

    def findBestPath(self, endState, fullTransitionGraph=None):
        """Find the best path from startState to endState.

        Args:
            endState (str): Ending state for finding the best path.
            fullTransitionGraph (networkx.DiGraph, optional): Pre-built transition graph.

        Returns:
            list: List of states representing the best path.
        """
        if fullTransitionGraph is None:
            fullTransitionGraph = self.makeFullTransitionGraph()
        bestPath = self.getBestPath(endState, fullTransitionGraph)
        return bestPath

    def makeFullTransitionGraph(self):
        """Create a full Transition graph based on stateTransitions.

        Returns:
            networkx.DiGraph: Full Transition graph.
        """
        transitionGraph = nx.DiGraph()
        for target, transitions in self.stateTransitions.items():
            for source in transitions.keys():
                forward = 0
                backward = 0
                if source in self.stateTransitions[target]:
                    forward = self.stateTransitions[target][source]
                if target in self.stateTransitions.get(source, {}):
                    backward = self.stateTransitions[source][target]
                totalWeight = forward - backward
                transitionGraph.add_edge(source, target, weight=totalWeight)
        return transitionGraph

    def getBestPath(self, endState, fullTransitionGraph):
        """Get the best path using Dijkstra's algorithm with a maximum depth.

        Args:
            endState (str): Ending state for finding the best path.
            fullTransitionGraph (networkx.DiGraph): Pre-built transition graph.

        Returns:
            list: List of states representing the best path.
        """
        heap = [((0, 0), [self.startState], float('inf'), 0)]
        bestPath = []
        bestScore = float('-inf')
        bestAverage = float('-inf')
        max_depth = self.maxDepthBest

        while heap:
            (priority, path, currentMinWeight, totalWeight) = heapq.heappop(heap)
            currentState = path[-1]
            current_depth = len(path) - 1

            if currentState == endState:
                currentAverage = totalWeight / len(path) if len(path) > 0 else 0
                if currentMinWeight > bestScore or (currentMinWeight == bestScore and currentAverage > bestAverage):
                    bestPath = path
                    bestScore = currentMinWeight
                    bestAverage = currentAverage
            elif current_depth == max_depth and bestPath == [] and max_depth < self.maxDepthPossible:
                max_depth = max_depth + 1
            
            if current_depth < max_depth:
                for neighbor in fullTransitionGraph.neighbors(currentState):
                    if neighbor not in path:
                        weight = fullTransitionGraph[currentState][neighbor]['weight']
                        new_min, new_total = self.calculateIncrementalScore(currentMinWeight, totalWeight, weight)
                        
                        if new_min < bestScore:
                            continue
                        
                        new_path = path + [neighbor]
                        new_average = new_total / len(new_path)
                        heapq.heappush(heap, ((-new_min, -new_average), new_path, new_min, new_total))
        return bestPath

    def calculateIncrementalScore(self, currentMinWeight, totalWeight, edgeWeight):
        """Calculate the next step minimum and total weights incrementally.

        Args:
            currentMinWeight (float): The current minimum weight along the path.
            totalWeight (float): The accumulated weight of the path.
            edgeWeight (float): The weight of the newly traversed edge.

        Returns:
            tuple: Updated (minWeight, totalWeight).
        """
        if currentMinWeight < 0 or edgeWeight < 0:
            return float('-inf'), totalWeight + edgeWeight
        return min(currentMinWeight, edgeWeight), totalWeight + edgeWeight

    def getPathScore(self, path, graph):
        """Calculate the minimum and mean weights for a given path.

        Args:
            path (list): List of states representing the path.
            graph (networkx.DiGraph): Transition graph.

        Returns:
            tuple: Tuple of (minWeight, meanWeight).
        """
        minWeight = float('inf')
        totalWeight = 0

        for i in range(len(path) - 1):
            weight = graph[path[i]][path[i + 1]]['weight']
            minWeight = min([minWeight, weight])
            totalWeight += weight

            if minWeight < 0:
                minWeight = float('-inf')
                return minWeight, 0

        meanWeight = totalWeight / len(path)
        return minWeight, meanWeight