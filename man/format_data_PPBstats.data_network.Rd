% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_data_PPBstats.data_network.R
\name{format_data_PPBstats.data_network}
\alias{format_data_PPBstats.data_network}
\title{Check and format the data to be used by PPBstats functions for network analyses}
\usage{
format_data_PPBstats.data_network(
  data,
  network_part = c("unipart", "bipart"),
  network_split = c("germplasm", "relation_year_start"),
  vertex_type = NULL
)
}
\arguments{
\item{data}{The data frame to format, see details.}

\item{network_part}{element of the network, it can be "unipart" or "bipart"}

\item{network_split}{For network_part = "unipart" and 
vertex_type = "location", split of the data that can be "germplasm" or "relation_year_start"}

\item{vertex_type}{\itemize{
\item for unipart network : "seed_lots" or "location"
\item for bipart network : c("germplasm", "location")
}}
}
\value{
It returns a igraph object coming from igraph::graph_from_data_frame().

For unipart network on seed lots, it a list of one element

For unipart network on location
\itemize{
 \item for network_split = "germplasm", 
  it returns a list with as many elements as germplam in the data
  as well as all germplasms merged in the first element of the list.
  \item for network_split = "relation_year_start", 
  it returns a list with as many elements as year in the data
  as well as all years merged in the first element of the list.
  }

For bipart network, it returns a list with as many elements as year in the data 
as well as all years merged in the first element of the list.
If no year are provided into the data, all information are merged.
}
\description{
\code{format_data_PPBstats} checks and formats the data to be used by PPBstats functions for network analyses
}
\details{
The data frame are different regarding type of network
 \itemize{
  \item for unipart network, two vertex_type are possible :
  \itemize{
   \item "seed_lots" : the data must have the following columns : 
   \itemize{
    \item "seed_lot_parent" : name of the seed lot parent in the relation
    \item "seed_lot_child" ; name of the seed lots child in the relation
    \item "relation_type" : the type of relation between the seed lots
    \item "relation_year_start" : the year when the relation starts
    \item "relation_year_end" : the year when the relation stops
    \item "germplasm_parent" : the germplasm associated to the seed lot father
    \item "location_parent" : the location associated to the seed lot father
    \item "year_parent" : represents the year of the last relation event of the seed lot father
    \item "germplasm_child" : the germplasm associated to the seed lot child
    \item "location_child" : the location associated to the seed lot child
    \item "year_child" : represents the year of the last relation event of the seed lot child
   }
   
   It can have in option : "alt_parent", "long_parent", "lat_parent",
   "alt_child", "long_child", "lat_child" to get map representation
   
   It can have supplementary variables with tags "_parent", "_child" or "_relation".
   
   \item "location" that represents each diffusion between location : the data can have two formats:
   \itemize{
    \item the same format than for unipart network and vertex_type = seed_lots
    \item the following columns (explained above): 
    "location_parent", "location_child"
    "relation_year_start", "relation_year_end"
    It can have in option : "germplasm_parent", "year_parent",
    "germplasm_child", "year_child"
    It can have in option : "alt_parent", "long_parent", "lat_parent",
    "alt_child", "long_child", "lat_child" to get map representation
    }
   }
  
  \item for bipartite network where a vertex can be a location or a germplasm, the data can have two formats:
  \itemize{
   \item the same format than for unipart network and vertex_type = seed_lots. 
   In this case, relation type diffusion or reproduction are kept.
   \item the following columns : "germplasm", "location", "year"
   It can have in option : "alt", "long", "lat" to get map representation
  }
 }

See the book for more details \href{https://priviere.github.io/PPBstats_book/introduction.html#data-network}{here}.
}
\seealso{
\code{\link{format_data_PPBstats}}
}
\author{
Pierre Riviere
}
