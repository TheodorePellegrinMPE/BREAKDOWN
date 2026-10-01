import networkx as nx
from graphviz import Digraph
import os
import concurrent.futures
from breakdown.transition_best_path import TransitionBestPath
from breakdown.auto_detect_utils import ensureFolderExists, getFragmentNames

class MoleculeTransitionGraph:
    def __init__(self, folder, stateNames, stateTransitions, stateTimes, netTransitions,
                 startState, numStatesToInclude, numFragStates, maxDepthBest, maxDepthPossible, states):
        """Initialize MoleculeTransitionGraph object.

        Args:
            folder (str): Folder name.
            stateTransitions (dict): Dictionary representing state transitions.
            stateTimes (dict): Dictionary mapping state keys to times.
            netTransitions (dict): Dictionary representing net transitions.
            startState (str): Name of start state isomer.
            numStatesToInclude (int): Number of top states to include in the flow chart.
            numFragStates (int): Number of scattered states to include in the flow chart.
            maxDepthBest (int): Depth to which the best path will be looked for.
            maxDepthPossible (int): Depth to which the any path will be looked for.
            states (list): List of states that must be visited.
        """
        self.saveFolder = os.path.abspath(folder)
        self.stateNames = stateNames
        self.stateTransitions = stateTransitions
        self.stateTimes = stateTimes
        self.netTransitions = netTransitions
        self.startState = startState
        self.numStatesToInclude = numStatesToInclude
        self.numFragStates = numFragStates
        self.manditoryStates = list(states)
        self.fullTransitionGraph = nx.DiGraph()
        self.bestPathObject = TransitionBestPath(self.stateTransitions, startState, maxDepthBest, maxDepthPossible)
        self.cSubstringCounts = self.precomputeCSubstrings()

    def precomputeCSubstrings(self):
        """Precompute the total count of 'C' substrings mapped directly to target values.

        Returns:
            dict: Mapping of state names to precomputed counts.
        """
        counts = {}
        for key, val in self.stateNames.items():
            count = sum(1 for substring in key.split('.') if 'C' in substring)
            counts[val] = counts.get(val, 0) + count
        return counts

    def makeTransitionGraph(self):
        """Create and save Transition graphs for top states using threading."""
        flowChartEndStates = self.getFlowChartEndStates()
        ensureFolderExists(f'{self.saveFolder}/transitionGraphs')
        maxWorkers = self.numStatesToInclude + self.numFragStates + len(self.manditoryStates)
        sharedGraph = self.bestPathObject.makeFullTransitionGraph()
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=maxWorkers) as executor:
            futures = []
            for endState in flowChartEndStates:
                futures.append(executor.submit(self.makePathGraphToEndState, endState, sharedGraph))

            concurrent.futures.wait(futures)
        for future in futures:
            future.result()
        fullDirectionalGraph = self.makeDiGraph(self.fullTransitionGraph)
        self.saveVisitedStates(self.fullTransitionGraph)
        fullDirectionalGraph.render(f'{self.saveFolder}/transitionGraphs/full_transitions_graph')

    def getFlowChartEndStates(self):
        """Get the most common states and the most common states where the molecule has blown up.

        Returns:
            list: The names of the states to be used to make the flow charts.
        """
        topStates = sorted(self.stateTimes,
                           key=lambda state: self.stateTimes[state],
                           reverse=True)[:self.numStatesToInclude]
        fragmentedStates = [state for state in self.stateTimes.keys() if (self.countCSubstrings(state) >= 4 and
                                                                          self.stateTimes[state] > 30)]
        numFragStates = min(self.numFragStates, len(fragmentedStates))
        flowChartEndStates = topStates + fragmentedStates[:numFragStates] + self.manditoryStates
        return flowChartEndStates
    
    def countCSubstrings(self, value):
        """Count the substrings containing the letter 'C' in the keys associated with a given value.

        Args:
            value (str): The value to search for in the dictionary.

        Returns:
            int: The count of substrings containing 'C' in the keys.
        """
        return self.cSubstringCounts.get(value, 0)

    def makePathGraphToEndState(self, endState, sharedGraph):
        """Create and save Transition graph for a specific end state.

        Args:
            endState (str): End state key.
            sharedGraph (networkx.DiGraph): Pre-built shared transition graph.
        """
        bestPath = self.bestPathObject.findBestPath(endState, sharedGraph)
        transitionGraph = nx.DiGraph()
        transitionGraph.add_nodes_from(bestPath)
        self.fullTransitionGraph.add_nodes_from(bestPath)
        for i in range(len(bestPath) - 1):
            self.addEdges(bestPath, transitionGraph, i)
        directionalGraph = self.makeDiGraph(transitionGraph)
        directionalGraph.render(f'{self.saveFolder}/transitionGraphs/transitions_graph_{endState}')

    def makeDiGraph(self, transitionGraph):
        """Create a directed graph using Graphviz.

        Args:
            transitionGraph (networkx.DiGraph): Transition graph.

        Returns:
            graphviz.Digraph: Directed graph.
        """
        directionalGraph = Digraph(format='pdf')
        directionalGraph.attr(rankdir='TB')
        for state in transitionGraph.nodes():
            directionalGraph = self.addStateToDiGraph(directionalGraph, state)
        for source, target, data in transitionGraph.edges(data=True):
            label = str(data['weight'])
            directionalGraph.edge(source, target, label=label)
        return directionalGraph

    def addStateToDiGraph(self, directionalGraph, state):
        """Add a state to the directed graph.

        Args:
            directionalGraph (graphviz.Digraph): Directed graph.
            state (str): State key.

        Returns:
            graphviz.Digraph: Updated directed graph.
        """
        time = self.stateTimes[state]
        label = f'{state}\n{time}'
        imagePath = f'{self.saveFolder}/drawings/{state}.png'
        if os.path.exists(imagePath):
            directionalGraph.node(state, label=label, image=imagePath, imagepos='tc', labelloc='b', height='5')
            directionalGraph.attr('node', shape='box', margin='0.1,0.1')
        else:
            directionalGraph.node(state, label=label, labelloc='b', height='5')
        return directionalGraph

    def saveVisitedStates(self, transitionGraph):
        """Save all nodes in the transition graph as a txt file.

        Args:
            transitionGraph (nx.DiGraph): DiGraph object containing all visited states.
        """
        isomerNamesSet = set()
        nodes = set(transitionGraph.nodes())
        for key, name in self.stateNames.items():
            if name in nodes:
                fragments, isomerNames = getFragmentNames(key, name)
                isomerNamesSet.update(isomerNames)
        with open(f'{self.saveFolder}/visited_states.txt', 'w') as file:
            for isomerName in isomerNamesSet:
                file.write(f"{isomerName}\n")

    def addEdges(self, bestPath, transitionGraph, i):
        """Add edges to the Transition graph.

        Args:
            bestPath (list): List of states representing the best path.
            transitionGraph (networkx.DiGraph): Transition graph.
            i (int): Index for bestPath.
        """
        fromState = bestPath[i]
        toState = bestPath[i + 1]
        forward = self.stateTransitions[toState][fromState]
        transitionGraph.add_edge(fromState, toState, weight=forward)
        self.fullTransitionGraph.add_edge(fromState, toState, weight=forward)
        if toState in self.stateTransitions.get(fromState, {}):
            backward = self.stateTransitions[fromState][toState]
            transitionGraph.add_edge(toState, fromState, weight=backward)
            self.fullTransitionGraph.add_edge(toState, fromState, weight=backward)
    
    def makeDirectTransitionGraph(self):
        """Creates a transition graph that includes the start state and all states directly
            connected to it.
        """
        transitionGraph = nx.DiGraph()
        transitionGraph.add_node(self.startState)
        for state in self.netTransitions:
            for fromState in self.netTransitions[state]:
                if fromState == self.startState:
                    forward = self.netTransitions[state][self.startState]
                    if forward != 0:
                        transitionGraph.add_node(state)
                        transitionGraph.add_edge(self.startState, state, weight=forward)
        ensureFolderExists(f'{self.saveFolder}/transitionGraphs')
        directionalGraph = Digraph(format='pdf')
        for state in transitionGraph.nodes():
            directionalGraph = self.addStateToDiGraph(directionalGraph, state)
        for edge in transitionGraph.edges(data=True):
            directionalGraph.edge(edge[0], edge[1], label=str(edge[2]['weight']))
        directionalGraph.render(f'{self.saveFolder}/transitionGraphs/direct_transitions_{self.startState}')
        self.saveVisitedStates(transitionGraph)