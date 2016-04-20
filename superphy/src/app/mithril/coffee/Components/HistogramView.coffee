class HistogramView
    constructor: () ->
        return @

    init: (searchResults) =>
    ## search results is an object
        # console.log("object", searchResults)
        # console.log("object keys", Object.keys(searchResults))
        genomes = Object.keys(searchResults)
        genes = Object.keys(searchResults[genomes[0]])

        view = @

        (el, isInitialized, ctx) ->
            #console.log(ctx)
            if not isInitialized
                view._create_histrogram(genomes, genes, searchResults, parentElem, elID)


    create_histogram: () =>
        # Create empty histogram
        margin = {top: 40, right: 30, bottom: 40, left: 30}
        @width = 300 - margin.left - margin.right
        @height = 250 - margin.top - margin.bottom

        bins = [
            {'val': 0, 'key': '0'}
            {'val': 1, 'key': '1'}
            {'val': 2, 'key': '2'}
            {'val': 3, 'key': '3'}
            {'val': 4, 'key': '4'}
            {'val': 5, 'key': '>=5'}
        ]

        @x = d3.scale.ordinal()
            .domain(bins.map((d) -> d.val))
            .rangeRoundBands([0, @width], .05)

        @x2 = d3.scale.ordinal()
            .domain(bins.map((d) -> d.key))
            .rangeRoundBands([0, @width], .05)

        @xAxis = d3.svg.axis()
            .scale(@x2)
            .orient("bottom")

        @histogram = d3.layout.histogram()
            .bins([0,1,2,3,4,5,6])

        @parentElem.append("<div id='#{@elID}' class='#{@cssClass}'></div>")
      
        @canvas = d3.select("##{@elID}").append("svg")
            .attr("width", @width + margin.left + margin.right)
            .attr("height", @height + margin.top + margin.bottom)
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
      
        @formatCount = d3.format(",.0f");
    
        @canvas.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + @height + ")")
            .call(@xAxis)
            .append("text")
            .attr("dy", ".75em")
            .attr("y", 23)
            .attr("x", @width / 2)
            .attr("text-anchor", "middle")
            .text( 'Number of Alleles')


