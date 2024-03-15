# NeGA

There are lots of things to do(will do it a long time later?).
The neighbor status is incorrectly set.
The probability of state update is determined by neural network. 
The general idea is not yet mature. 
Just pseudocode below~

```
// Initializing the Metacellular Automata Model.
CREATE ARRAY cells [width][height]  
FOR i FROM 0 TO width DO  
    FOR j FROM 0 TO height DO  
            cells[i][j] = initialize_cell_state(i, j) 
	    // Metacells are initialized according to the specific problem.  
    END FOR  
END FOR  
  
// Main Loop.
WHILE True DO  // To simplify the pseudo-code, 
               // we actually use time span and step size.
    FOR i FROM 0 TO width DO  
        FOR j FROM 0 TO height DO  
            // Update the state according to the rules of the metacells.
            IF cells[i][j] == 0  
		// If the current metacell is inorganic environment.
                neighbor_state = get_random_neighbor_state(cells[i][j])
		// Randomly read the state of one metacell around the current metacell.
                cells[i][j] = update_cell_state(neighbor_state, cells[i][j])
		// Based on the state of its neighbors and its own state, 
// as well as the effects of random events, it finally decides its own state
// at the next moment in a comprehensive manner.
                IF cells[i][j] != 0
                    cells[i][j] = disaster(stochastic_disaster_cycle,cells[i][j],
                                                        disaster_resilience[i][j])
                END IF
                // Considering the effect of stochastic disaster cycles on metacell.
           
	    ... // Partial omission.

            ELSEIF cells[i][j] == 5
                // If the current metacell is toppredator.
                neighbor_state = get_random_neighbor_state(cells[i][j])
                cells[i][j] = update_cell_state(neighbor_state, cells[i][j])
                IF cells[i][j] != 0
                    cells[i][j] = disaster(stochastic_disaster_cycle,cells[i][j],
                                                        disaster_resilience[i][j])
                END IF
            END IF  
        END FOR  
    END FOR  
END WHILE
```
