class Singleton
    model: () =>
    constructor: () ->
        @model()
    @get: ->
        @instance ?= new @()
        return @instance
    @getView: () ->
        return @get().view()
    view: () ->    

class Page_Template extends Singleton
    model: () ->
    controller: () ->
    view: () ->
        return [
            Header.getView()
        ]