class Singleton
    model: () ->
    constructor: () ->
        @model()
    @get: ->
        @instance ?= new @()
        return @instance
    @getView: () ->
        return @get().view()
    view: () ->

class PageTemplate extends Singleton
    model: () ->
    controller: () ->
    view: () ->
        return [
            Header.getView()
        ]