# secondary-structure-parser
Parser for ss3 and ss8 secondary structure predictions.
Files are parsed and validated. Statistics like secondary-structure distribution are calculated on-demand.
## Usage

### SecondaryStructureParser
```

from SecondaryStructureParser import SecondaryStructureParser
parser = SecondaryStructureParser('examples/1pazA.ss3')
print(parser.calculate_statistics())

>>> {'C': 0.4959349593495935, 'E': 0.34959349593495936, 'H': 0.15447154471544716}
```

or

`python SecondaryStructureParser example.ss3`

Individual entries can be accesed via `parser.parsed[residueId]`
### BatchParser

```
from BatchParser import BatchParser
batch_parser = BatchParser('examples')
print(batch_parser.calculate_statistics())
```

or

`python BatchParser example.ss3`

```
{'total number of results': 1, 
 'predictions': {'E': 0.34959349593495936, 'H': 0.15447154471544716, 'C': 0.4959349593495935}, 
 'summary': {'E': {'max': 0.34959349593495936, 'median': 0.34959349593495936, 'average': 0.34959349593495936, 'stdev': 0.0, 'min': 0.34959349593495936}, 
             'H': {'max': 0.15447154471544716, 'median': 0.15447154471544716, 'average': 0.15447154471544716, 'stdev': 0.0, 'min': 0.15447154471544716}, 
             'C': {'max': 0.4959349593495935, 'median': 0.4959349593495935, 'average': 0.4959349593495935, 'stdev': 0.0, 'min': 0.4959349593495935}}, 
 'examples/1pazA.ss3': {'E': 0.34959349593495936, 'H': 0.15447154471544716, 'C': 0.4959349593495935}}
 ```  

Individual files can be accesed via `batch_parser.parsed[filename]`     
     