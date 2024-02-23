
#' @noRd
#' @importFrom ggplot2 guide_train
guide_train.prism_axis <- function(guide, scale, aesthetic = NULL) {

  aesthetic <- aesthetic %||% scale$aesthetics[1]
  breaks <- scale$get_breaks()

  empty_ticks <- base::data.frame(
    aesthetic = numeric(0),
    .value    = numeric(0),
    .label    = character(0)
  )
  names(empty_ticks) <- c(aesthetic, ".value", ".label")

  if (length(intersect(scale$aesthetics, guide$available_aes)) == 0) {
    warn(glue(
      "axis guide needs appropriate scales: ",
      glue_collapse(guide$available_aes, ", ", last = " or ")
    ))
    guide$key <- empty_ticks
  } else if (length(breaks) == 0) {
    guide$key <- empty_ticks
  } else {
    mapped_breaks <- if (scale$is_discrete()) scale$map(breaks) else breaks
    ticks <- base::data.frame(setNames(list(mapped_breaks), aesthetic))
    ticks$.value <- breaks
    ticks$.label <- scale$get_labels(breaks)

    guide$key <- ticks[is.finite(ticks[[aesthetic]]), ]
  }

  guide$name <- paste0(guide$name, "_", aesthetic)
  guide$hash <- digest(list(guide$title, guide$key$.value,
                            guide$key$.label, guide$name, guide$is_major))
  guide
}

#' @noRd
#' @importFrom ggplot2 guide_transform
guide_transform.prism_axis <- function(guide, coord, panel_params) {

  if (is.null(guide$position) || nrow(guide$key) == 0) {
    return(guide)
  }

  aesthetics <- names(guide$key)[!grepl("^\\.", names(guide$key))]

  if (all(c("x", "y") %in% aesthetics)) {
    guide$key <- coord$transform(guide$key, panel_params)
  } else {
    other_aesthetic <- setdiff(c("x", "y"), aesthetics)
    override_value  <- if(guide$position %in% c("bottom", "left")) -Inf else Inf
    guide$key[[other_aesthetic]] <- override_value

    guide$key <- coord$transform(guide$key, panel_params)

    warn_for_guide_position(guide)
  }
  guide
}

warn_for_guide_position <- function(guide) {
  # This is trying to catch when a user specifies a position perpendicular
  # to the direction of the axis (e.g., a "y" axis on "top").
  # The strategy is to check that two or more unique breaks are mapped
  # to the same value along the axis.
  breaks_are_unique <- !duplicated(guide$key$.value)
  empty <- is.null(guide$key) || prod(dim(guide$key)) == 0 ||
    inherits(guide$key, "waiver")
  if (empty || sum(breaks_are_unique) == 1) {
    return()
  }

  if (guide$position %in% c("top", "bottom")) {
    position_aes <- "x"
  } else if(guide$position %in% c("left", "right")) {
    position_aes <- "y"
  } else {
    return()
  }

  if (length(unique(guide$key[[position_aes]][breaks_are_unique])) == 1) {
    warn(c(
      "Position guide is perpendicular to the intended axis",
      "i" = "Did you mean to specify a different guide `position`?"
    ))
  }
}

#' @noRd
#' @importFrom ggplot2 guide_geom
guide_geom.prism_axis <- function(guide, layers, ...) {
  guide
}

#' @noRd
#' @importFrom ggplot2 guide_merge
guide_merge.prism_axis <- function(guide, new_guide) {
  if (!inherits(new_guide, c("guide_none", "GuideNone"))) {
    warn(c(
      "Discarding guide on merge",
      "i" = "Do you have more than one guide with the same position?"
    ))
  }
}

#' Axis guide with brackets
#'
#' This guide turns the axis into brackets drawn around each axis label.
#'
#' The number of brackets can be adjusted using the `breaks`
#' argument in `scale_(x|y)_continuous()` or `scale_(x|y)_discrete()`.
#'
#' @inheritParams ggplot2::guide_axis
#' @param width `numeric`. Controls the width of the bracket. Try
#' values between 0 and 1.
#' @param outside `logical`. Default is `TRUE` and brackets point
#' outwards. If `FALSE` the bracket crossbar is moved so the ticks appear
#' to point inwards towards the plotting area.
#'
#' @return Returns a \emph{prism_bracket} guide class object.
#'
#' @example inst/examples/ex-guide_prism_bracket.R
#'
#' @noRd
guide_prism_bracket <- function(title = waiver(), check.overlap = FALSE,
                                angle = NULL, n.dodge = 1, order = 0,
                                position = waiver(), width = NULL,
                                outside = TRUE) {
  if (packageVersion("ggplot2") < "3.3.0") {
    stop("ggplot2 >= 3.3.0 needed for this function.", call. = FALSE)
  }

  structure(
    list(
      title = title,

      # customizations
      check.overlap = check.overlap,
      angle = angle,
      n.dodge = n.dodge,

      # general
      order = order,
      position = position,

      # parameter
      available_aes = c("x", "y"),

      # custom
      width = width,
      outside = outside,

      name = "axis"
    ),
    class = c("guide", "prism_bracket", "prism_axis", "axis")
  )
}

#' @keywords internal
#' @noRd
guide_gengrob.prism_bracket <- function(guide, theme) {
  aesthetic <- names(guide$key)[!grepl("^\\.", names(guide$key))][1]

  draw_prism_bracket(
    break_positions = guide$key[[aesthetic]],
    break_labels = guide$key$.label,
    axis_position = guide$position,
    theme = theme,
    check.overlap = guide$check.overlap,
    angle = guide$angle,
    n.dodge = guide$n.dodge,
    width = guide$width,
    outside = guide$outside
  )
}


#' Grob for bracket axes
#'
#' @description Grob for bracket axes.
#'
#' @param break_positions Position of bracket center and labels
#' @param break_labels Labels between ticks
#' @param axis_position Position of axis (top, bottom, left or right)
#' @param theme A complete \code{\link[ggplot2]{theme}} object
#' @param check.overlap Silently remove overlapping labels,
#'   (recursively) prioritizing the first, last, and middle labels.
#' @param angle Compared to setting the angle in
#'   \code{\link[ggplot2]{theme}} / `element_text`,
#'   this also uses some heuristics to automatically pick the `hjust` and
#'   `vjust` that you probably want.
#' @param n.dodge The number of rows (for vertical axes) or columns (for
#'   horizontal axes) that should be used to render the labels. This is
#'   useful for displaying labels that would otherwise overlap.
#' @param width `numeric`. Controls the width of the bracket. Try
#' values between 0 and 1.
#' @param outside `logical`. Default is `TRUE` and brackets point
#' outwards. If `FALSE` the bracket crossbar is moved so the ticks appear
#' to point inwards towards the plotting area.
#' @keywords internal
draw_prism_bracket <- function(break_positions, break_labels, axis_position,
                               theme, check.overlap = FALSE, angle = NULL,
                               n.dodge = 1, width = NULL,
                               outside = TRUE) {

  axis_position <- match.arg(axis_position, c("top", "bottom", "right", "left"))
  aesthetic <- if (axis_position %in% c("top", "bottom")) "x" else "y"

  # resolve elements
  line_element_name <- paste0("axis.line.", aesthetic, ".", axis_position)
  tick_element_name <- paste0("axis.ticks.", aesthetic, ".", axis_position)
  tick_length_element_name <- paste0("axis.ticks.length.", aesthetic, ".", axis_position)
  label_element_name <- paste0("axis.text.", aesthetic, ".", axis_position)

  line_element <- calc_element(line_element_name, theme)
  tick_element <- calc_element(tick_element_name, theme)
  tick_length <- calc_element(tick_length_element_name, theme)
  label_element <- calc_element(label_element_name, theme)

  # override label element parameters for rotation
  if (inherits(label_element, "element_text")) {
    label_overrides <- axis_label_element_overrides(axis_position, angle)
    # label_overrides is an element_text, but label_element may not be;
    # to merge the two elements, we just copy angle, hjust, and vjust
    # unless their values are NULL
    if (!is.null(label_overrides$angle)) {
      label_element$angle <- label_overrides$angle
    }
    if (!is.null(label_overrides$hjust)) {
      label_element$hjust <- label_overrides$hjust
    }
    if (!is.null(label_overrides$vjust)) {
      label_element$vjust <- label_overrides$vjust
    }
  }

  # conditionally set parameters that depend on axis orientation
  is_vertical <- axis_position %in% c("left",  "right")

  position_dim <- if (is_vertical) "y" else "x"
  non_position_dim <- if (is_vertical) "x" else "y"
  position_size <- if (is_vertical) "height" else "width"
  non_position_size <- if (is_vertical) "width" else "height"
  gtable_element <- if (is_vertical) gtable_row else gtable_col
  measure_gtable <- if (is_vertical) gtable_width else gtable_height
  measure_labels_non_pos <- if (is_vertical) grobWidth else grobHeight

  # conditionally set parameters that depend on which side of the panel
  # the axis is on
  is_second <- axis_position %in% c("right", "top")

  tick_direction <- if (is_second) 1 else -1
  non_position_panel <- if (is_second) unit(0, "npc") else unit(1, "npc")
  tick_coordinate_order <- if (is_second) c(2, 1) else c(1, 2)

  # conditionally set the gtable ordering
  labels_first_gtable <- axis_position %in% c("left", "top") # position in gtable

  # set common parameters
  n_breaks <- length(break_positions)
  opposite_positions <- c("top" = "bottom", "bottom" = "top",
                          "right" = "left", "left" = "right")
  axis_position_opposite <- unname(opposite_positions[axis_position])

  # autocalculate bracket width based on number of breaks if missing
  # best for discrete axes, bad for continuous axes
  if (is.null(width)) {
    if (n_breaks == 1) {
      width <- 0.75
    }
    else if (n_breaks > 1) {
      width <- (0.8 + 0.01 * n_breaks) / n_breaks
    }
  }

  # draw elements
  half_bracket <- width / 2

  lines_grob <- exec(
    element_grob, line_element,
    !!position_dim := unit.c(
      unit(
        sort(c(break_positions - half_bracket,
               break_positions + half_bracket)), "native"
      )
    ),
    !!non_position_dim := if (outside) {
      rep(
        unit.c(non_position_panel, non_position_panel),
        times = n_breaks
      )
    } else {
      rep(
        unit.c(non_position_panel + (tick_direction * tick_length),
               non_position_panel + (tick_direction * tick_length)),
        times = n_breaks
      )
    },
    id.lengths = rep(2, times = n_breaks)
  )

  if (n_breaks == 0) {
    return(
      absoluteGrob(
        gList(lines_grob),
        width = grobWidth(lines_grob),
        height = grobHeight(lines_grob)
      )
    )
  }

  # break_labels can be a list() of language objects
  if (is.list(break_labels)) {
    if (any(vapply(break_labels, is.language, logical(1)))) {
      break_labels <- do.call(expression, break_labels)
    } else {
      break_labels <- unlist(break_labels)
    }
  }

  # calculate multiple rows/columns of labels (which is usually 1)
  dodge_pos <- rep(seq_len(n.dodge), length.out = n_breaks)
  dodge_indices <- split(seq_len(n_breaks), dodge_pos)

  label_grobs <- lapply(dodge_indices, function(indices) {
    draw_axis_labels(
      break_positions = break_positions[indices],
      break_labels = break_labels[indices],
      label_element = label_element,
      is_vertical = is_vertical,
      check.overlap = check.overlap
    )
  })

  ticks_grob <- exec(
    element_grob, tick_element,
    !!position_dim := unit.c(
      rep(unit(break_positions - half_bracket, "native"), each = 2),
      rep(unit(break_positions + half_bracket, "native"), each = 2)
    ),
    !!non_position_dim := rep(
      unit.c(non_position_panel + (tick_direction * tick_length),
             non_position_panel)[tick_coordinate_order],
      times = n_breaks * 2
    ),
    id.lengths = rep(2, times = n_breaks * 2)
  )

  # create gtable
  non_position_sizes <- paste0(non_position_size, "s")
  label_dims <- do.call(unit.c, lapply(label_grobs, measure_labels_non_pos))
  grobs <- c(list(ticks_grob), label_grobs)
  grob_dims <- unit.c(tick_length, label_dims)

  if (labels_first_gtable) {
    grobs <- rev(grobs)
    grob_dims <- rev(grob_dims)
  }

  gt <- exec(
    gtable_element,
    name = "axis",
    grobs = grobs,
    !!non_position_sizes := grob_dims,
    !!position_size := unit(1, "npc")
  )

  # create viewport
  justvp <- exec(
    viewport,
    !!non_position_dim := non_position_panel,
    !!non_position_size := measure_gtable(gt),
    just = axis_position_opposite
  )

  absoluteGrob(
    gList(lines_grob, gt),
    width = gtable_width(gt),
    height = gtable_height(gt),
    vp = justvp
  )
}


#' Axis guide with minor ticks
#'
#' This guide is like the standard \code{\link[ggplot2]{guide_axis}}, but
#' with minor ticks.
#'
#' The number of minor ticks can be changed using the `minor_breaks`
#' argument. Control the length of minor ticks by setting
#' `prism.ticks.length` to a \code{\link[grid]{unit}} object using
#' \code{\link[ggplot2]{theme}}, for example:
#' `prism.ticks.length = unit(2, "pt")`. The major tick lengths
#' are adjusted using the standard `axis.ticks.length`.
#'
#' @inheritParams ggplot2::guide_axis
#'
#' @return Returns a \emph{prism_minor} guide class object.
#'
#' @example inst/examples/ex-guide_prism_minor.R
#'
#' @noRd
guide_prism_minor <- function(title = waiver(), check.overlap = FALSE,
                              angle = NULL, n.dodge = 1, order = 0,
                              position = waiver()) {
  if (packageVersion("ggplot2") < "3.3.0") {
    stop("ggplot2 >= 3.3.0 needed for this function.", call. = FALSE)
  }

  structure(
    list(
      title = title,

      # customizations
      check.overlap = check.overlap,
      angle = angle,
      n.dodge = n.dodge,

      # general
      order = order,
      position = position,

      # parameter
      available_aes = c("x", "y"),

      name = "axis"
    ),
    class = c("guide", "prism_minor", "prism_axis", "axis")
  )
}

#' @keywords internal
#' @noRd
guide_train.prism_minor <- function(guide, scale, aesthetic = NULL) {

  aesthetic <- aesthetic %||% scale$aesthetics[1]

  # define major and minor breaks
  major_breaks <- scale$get_breaks()
  major_breaks <- major_breaks[!is.na(major_breaks)]

  if (is.null(scale$minor_breaks)) {
    stop("No minor breaks exist, guide_prism_minor needs minor breaks to work")
  }
  minor_breaks <- setdiff(scale$get_breaks_minor(), major_breaks)

  # define all breaks
  breaks <- union(major_breaks, minor_breaks)

  # indicate which breaks are major
  is_major <- breaks %in% major_breaks

  empty_ticks <- base::data.frame(
    list(aesthetic = numeric(0), .value = numeric(0), .label = character(0))
  )
  names(empty_ticks) <- c(aesthetic, ".value", ".label")

  if (length(intersect(scale$aesthetics, guide$available_aes)) == 0) {
    warn(glue(
      "axis guide needs appropriate scales: ",
      glue_collapse(guide$available_aes, ", ", last = " or ")
    ))
    guide$key <- empty_ticks
  } else if (length(breaks) == 0) {
    guide$key <- empty_ticks
  } else {
    mapped_breaks <- if (scale$is_discrete()) scale$map(breaks) else breaks

    ticks <- base::data.frame(setNames(list(mapped_breaks), aesthetic))
    ticks$.value <- breaks
    # get major break labels and make minor breaks blank
    ticks$.label <- c(scale$get_labels(major_breaks),
                      rep("", times = length(minor_breaks)))

    guide$key <- ticks[is.finite(ticks[[aesthetic]]), ]
  }

  guide$name <- paste0(guide$name, "_", aesthetic)
  guide$is_major <- is_major
  guide$hash <- digest(list(guide$title, guide$key$.value,
                            guide$key$.label, guide$name, guide$is_major))
  guide
}

#' @keywords internal
#' @noRd
guide_gengrob.prism_minor <- function(guide, theme) {
  aesthetic <- names(guide$key)[!grepl("^\\.", names(guide$key))][1]

  draw_prism_minor(
    break_positions = guide$key[[aesthetic]],
    break_labels = guide$key$.label[guide$is_major],
    breaks_major = guide$is_major,
    axis_position = guide$position,
    theme = theme,
    check.overlap = guide$check.overlap,
    angle = guide$angle,
    n.dodge = guide$n.dodge
  )
}

#' Grob for axes with minor ticks
#'
#' @param break_positions position of ticks
#' @param break_labels labels at ticks
#' @param breaks_major logical vector indicating major ticks versus minor ticks
#' @param axis_position position of axis (top, bottom, left or right)
#' @param theme A complete \code{\link[ggplot2]{theme}} object
#' @param check.overlap silently remove overlapping labels,
#'   (recursively) prioritizing the first, last, and middle labels.
#' @param angle Compared to setting the angle in
#'   \code{\link[ggplot2]{theme}} / `element_text`,
#'   this also uses some heuristics to automatically pick the `hjust` and
#'   `vjust` that you probably want.
#' @param n.dodge The number of rows (for vertical axes) or columns (for
#'   horizontal axes) that should be used to render the labels. This is
#'   useful for displaying labels that would otherwise overlap.
#' @keywords internal
draw_prism_minor <- function(break_positions, break_labels, breaks_major,
                             axis_position, theme,
                             check.overlap = FALSE, angle = NULL, n.dodge = 1) {

  axis_position <- match.arg(axis_position, c("top", "bottom", "right", "left"))
  aesthetic <- if (axis_position %in% c("top", "bottom")) "x" else "y"

  # resolve elements
  line_element_name <- paste0("axis.line.", aesthetic, ".", axis_position)
  tick_element_name <- paste0("axis.ticks.", aesthetic, ".", axis_position)
  tick_length_element_name <- paste0("axis.ticks.length.", aesthetic, ".", axis_position)
  label_element_name <- paste0("axis.text.", aesthetic, ".", axis_position)
  prism_tick_element_name <- paste0("prism.ticks.", aesthetic, ".", axis_position)
  prism_tick_length_element_name <- paste0("prism.ticks.length.", aesthetic, ".", axis_position)

  line_element <- calc_element(line_element_name, theme)
  tick_element <- calc_element(tick_element_name, theme)
  tick_length <- calc_element(tick_length_element_name, theme)
  label_element <- calc_element(label_element_name, theme)
  prism_tick_element <- calc_element(prism_tick_element_name, theme)
  prism_tick_length <- calc_element(prism_tick_length_element_name, theme)

  # override label element parameters for rotation
  if (inherits(label_element, "element_text")) {
    label_overrides <- axis_label_element_overrides(axis_position, angle)
    # label_overrides is an element_text, but label_element may not be;
    # to merge the two elements, we just copy angle, hjust, and vjust
    # unless their values are NULL
    if (!is.null(label_overrides$angle)) {
      label_element$angle <- label_overrides$angle
    }
    if (!is.null(label_overrides$hjust)) {
      label_element$hjust <- label_overrides$hjust
    }
    if (!is.null(label_overrides$vjust)) {
      label_element$vjust <- label_overrides$vjust
    }
  }

  # conditionally set parameters that depend on axis orientation
  is_vertical <- axis_position %in% c("left",  "right")

  position_dim <- if (is_vertical) "y" else "x"
  non_position_dim <- if (is_vertical) "x" else "y"
  position_size <- if (is_vertical) "height" else "width"
  non_position_size <- if (is_vertical) "width" else "height"
  gtable_element <- if (is_vertical) gtable_row else gtable_col
  measure_gtable <- if (is_vertical) gtable_width else gtable_height
  measure_labels_non_pos <- if (is_vertical) grobWidth else grobHeight

  # conditionally set parameters that depend on which side of the panel
  # the axis is on
  is_second <- axis_position %in% c("right", "top")

  tick_direction <- if (is_second) 1 else -1
  non_position_panel <- if (is_second) unit(0, "npc") else unit(1, "npc")
  tick_coordinate_order <- if (is_second) c(2, 1) else c(1, 2)

  # conditionally set the gtable ordering
  labels_first_gtable <- axis_position %in% c("left", "top") # refers to position in gtable

  # set common parameters
  n_breaks <- length(break_positions[breaks_major])
  n_minor_breaks <- length(break_positions[!breaks_major])
  opposite_positions <- c("top" = "bottom", "bottom" = "top", "right" = "left", "left" = "right")
  axis_position_opposite <- unname(opposite_positions[axis_position])

  # draw elements
  line_grob <- exec(
    element_grob, line_element,
    !!position_dim := unit(c(0, 1), "npc"),
    !!non_position_dim := unit.c(non_position_panel, non_position_panel)
  )

  if (n_breaks == 0) {
    return(
      absoluteGrob(
        gList(line_grob),
        width = grobWidth(line_grob),
        height = grobHeight(line_grob)
      )
    )
  }

  # break_labels can be a list() of language objects
  if (is.list(break_labels)) {
    if (any(vapply(break_labels, is.language, logical(1)))) {
      break_labels <- do.call(expression, break_labels)
    } else {
      break_labels <- unlist(break_labels)
    }
  }

  # calculate multiple rows/columns of labels (which is usually 1)
  dodge_pos <- rep(seq_len(n.dodge), length.out = n_breaks)
  dodge_indices <- split(seq_len(n_breaks), dodge_pos)

  label_grobs <- lapply(dodge_indices, function(indices) {
    draw_axis_labels(
      break_positions = break_positions[breaks_major][indices],
      break_labels = break_labels[indices],
      label_element = label_element,
      is_vertical = is_vertical,
      check.overlap = check.overlap
    )
  })

  ticks_grob <- exec(
    element_grob, tick_element,
    !!position_dim := rep(unit(break_positions, "native"), each = 2),
    !!non_position_dim := unit.c(
      rep(
        unit.c(non_position_panel + (tick_direction * tick_length), non_position_panel)[tick_coordinate_order],
        times = n_breaks
      ),
      rep(
        # let minor ticks have a different length from major ticks
        unit.c(non_position_panel + (tick_direction * prism_tick_length), non_position_panel)[tick_coordinate_order],
        times = n_minor_breaks
      )
    ),
    id.lengths = rep(2, times = length(break_positions))
  )

  # create gtable
  non_position_sizes <- paste0(non_position_size, "s")
  label_dims <- do.call(unit.c, lapply(label_grobs, measure_labels_non_pos))
  grobs <- c(list(ticks_grob), label_grobs)
  grob_dims <- unit.c(tick_length, label_dims)

  if (labels_first_gtable) {
    grobs <- rev(grobs)
    grob_dims <- rev(grob_dims)
  }

  gt <- exec(
    gtable_element,
    name = "axis",
    grobs = grobs,
    !!non_position_sizes := grob_dims,
    !!position_size := unit(1, "npc")
  )

  # create viewport
  justvp <- exec(
    viewport,
    !!non_position_dim := non_position_panel,
    !!non_position_size := measure_gtable(gt),
    just = axis_position_opposite
  )

  absoluteGrob(
    gList(line_grob, gt),
    width = gtable_width(gt),
    height = gtable_height(gt),
    vp = justvp
  )
}

#' Offset axis guide
#'
#' This guide draws the axis only as wide as the outermost tick marks,
#' similar to offset axes from Prism.
#'
#' Control the length of the axis by adjusting the `breaks` argument in
#' `scale_(x|y)_continuous()` or `scale_(x|y)_discrete()`.
#'
#' @inheritParams ggplot2::guide_axis
#'
#' @return Returns a \emph{prism_offset} guide class object.
#'
#' @example inst/examples/ex-guide_prism_offset.R
#'
#' @export
guide_prism_offset <- function(title = waiver(), check.overlap = FALSE,
                               angle = NULL, n.dodge = 1, order = 0,
                               position = waiver()) {
  if (packageVersion("ggplot2") < "3.3.0") {
    stop("ggplot2 >= 3.3.0 needed for this function.", call. = FALSE)
  }

  structure(
    list(
      title = title,

      # customizations
      check.overlap = check.overlap,
      angle = angle,
      n.dodge = n.dodge,

      # general
      order = order,
      position = position,

      # parameter
      available_aes = c("x", "y"),

      name = "axis"
    ),
    class = c("guide", "prism_offset", "prism_axis", "axis")
  )
}

#' @keywords internal
#' @export
guide_gengrob.prism_offset <- function(guide, theme) {
  aesthetic <- names(guide$key)[!grepl("^\\.", names(guide$key))][1]

  draw_prism_offset(
    break_positions = guide$key[[aesthetic]],
    break_labels = guide$key$.label,
    axis_position = guide$position,
    theme = theme,
    check.overlap = guide$check.overlap,
    angle = guide$angle,
    n.dodge = guide$n.dodge
  )
}


#' Grob for offset axes
#'
#' @param break_positions position of ticks
#' @param break_labels labels at ticks
#' @param axis_position position of axis (top, bottom, left or right)
#' @param theme A complete \code{\link[ggplot2]{theme}} object
#' @param check.overlap silently remove overlapping labels,
#'   (recursively) prioritizing the first, last, and middle labels.
#' @param angle Compared to setting the angle in
#'   \code{\link[ggplot2]{theme}} / `element_text`,
#'   this also uses some heuristics to automatically pick the `hjust` and
#'   `vjust` that you probably want.
#' @param n.dodge The number of rows (for vertical axes) or columns (for
#'   horizontal axes) that should be used to render the labels. This is
#'   useful for displaying labels that would otherwise overlap.
#' @keywords internal
draw_prism_offset <- function(break_positions, break_labels, axis_position, theme,
                              check.overlap = FALSE, angle = NULL, n.dodge = 1) {

  axis_position <- match.arg(axis_position, c("top", "bottom", "right", "left"))
  aesthetic <- if (axis_position %in% c("top", "bottom")) "x" else "y"

  # resolve elements
  line_element_name <- paste0("axis.line.", aesthetic, ".", axis_position)
  tick_element_name <- paste0("axis.ticks.", aesthetic, ".", axis_position)
  tick_length_element_name <- paste0("axis.ticks.length.", aesthetic, ".", axis_position)
  label_element_name <- paste0("axis.text.", aesthetic, ".", axis_position)

  line_element <- calc_element(line_element_name, theme)
  tick_element <- calc_element(tick_element_name, theme)
  tick_length <- calc_element(tick_length_element_name, theme)
  label_element <- calc_element(label_element_name, theme)

  # override label element parameters for rotation
  if (inherits(label_element, "element_text")) {
    label_overrides <- axis_label_element_overrides(axis_position, angle)
    # label_overrides is an element_text, but label_element may not be;
    # to merge the two elements, we just copy angle, hjust, and vjust
    # unless their values are NULL
    if (!is.null(label_overrides$angle)) {
      label_element$angle <- label_overrides$angle
    }
    if (!is.null(label_overrides$hjust)) {
      label_element$hjust <- label_overrides$hjust
    }
    if (!is.null(label_overrides$vjust)) {
      label_element$vjust <- label_overrides$vjust
    }
  }

  # conditionally set parameters that depend on axis orientation
  is_vertical <- axis_position %in% c("left",  "right")

  position_dim <- if (is_vertical) "y" else "x"
  non_position_dim <- if (is_vertical) "x" else "y"
  position_size <- if (is_vertical) "height" else "width"
  non_position_size <- if (is_vertical) "width" else "height"
  gtable_element <- if (is_vertical) gtable_row else gtable_col
  measure_gtable <- if (is_vertical) gtable_width else gtable_height
  measure_labels_non_pos <- if (is_vertical) grobWidth else grobHeight

  # conditionally set parameters that depend on which side of the panel
  # the axis is on
  is_second <- axis_position %in% c("right", "top")

  tick_direction <- if (is_second) 1 else -1
  non_position_panel <- if (is_second) unit(0, "npc") else unit(1, "npc")
  tick_coordinate_order <- if (is_second) c(2, 1) else c(1, 2)

  # conditionally set the gtable ordering
  labels_first_gtable <- axis_position %in% c("left", "top") # refers to position in gtable

  # set common parameters
  n_breaks <- length(break_positions)
  opposite_positions <- c("top" = "bottom", "bottom" = "top", "right" = "left", "left" = "right")
  axis_position_opposite <- unname(opposite_positions[axis_position])

  # draw elements
  line_grob <- exec(
    element_grob, line_element,
    !!position_dim := unit(c(min(break_positions),
                             max(break_positions)), "npc"),
    !!non_position_dim := unit.c(non_position_panel, non_position_panel)
  )

  if (n_breaks == 0) {
    return(
      absoluteGrob(
        gList(line_grob),
        width = grobWidth(line_grob),
        height = grobHeight(line_grob)
      )
    )
  }

  # break_labels can be a list() of language objects
  if (is.list(break_labels)) {
    if (any(vapply(break_labels, is.language, logical(1)))) {
      break_labels <- do.call(expression, break_labels)
    } else {
      break_labels <- unlist(break_labels)
    }
  }

  # calculate multiple rows/columns of labels (which is usually 1)
  dodge_pos <- rep(seq_len(n.dodge), length.out = n_breaks)
  dodge_indices <- split(seq_len(n_breaks), dodge_pos)

  label_grobs <- lapply(dodge_indices, function(indices) {
    draw_axis_labels(
      break_positions = break_positions[indices],
      break_labels = break_labels[indices],
      label_element = label_element,
      is_vertical = is_vertical,
      check.overlap = check.overlap
    )
  })

  ticks_grob <- exec(
    element_grob, tick_element,
    !!position_dim := rep(unit(break_positions, "native"), each = 2),
    !!non_position_dim := rep(
      unit.c(non_position_panel + (tick_direction * tick_length), non_position_panel)[tick_coordinate_order],
      times = n_breaks
    ),
    id.lengths = rep(2, times = n_breaks)
  )

  # create gtable
  non_position_sizes <- paste0(non_position_size, "s")
  label_dims <- do.call(unit.c, lapply(label_grobs, measure_labels_non_pos))
  grobs <- c(list(ticks_grob), label_grobs)
  grob_dims <- unit.c(tick_length, label_dims)

  if (labels_first_gtable) {
    grobs <- rev(grobs)
    grob_dims <- rev(grob_dims)
  }

  gt <- exec(
    gtable_element,
    name = "axis",
    grobs = grobs,
    !!non_position_sizes := grob_dims,
    !!position_size := unit(1, "npc")
  )

  # create viewport
  justvp <- exec(
    viewport,
    !!non_position_dim := non_position_panel,
    !!non_position_size := measure_gtable(gt),
    just = axis_position_opposite
  )

  absoluteGrob(
    gList(line_grob, gt),
    width = gtable_width(gt),
    height = gtable_height(gt),
    vp = justvp
  )
}

#' Offset axis guide with minor ticks
#'
#' This guide draws the axis only as wide as the outermost tick marks,
#' similar to offset axes from Prism. It also adds minor ticks.
#'
#' Control the length of the axis by adjusting the `breaks` argument in
#' `scale_(x|y)_continuous()` or `scale_(x|y)_discrete()`. Similarly,
#' the number of minor ticks can be changed using the `minor_breaks`
#' argument.
#'
#' Control the length of minor ticks by setting `prism.ticks.length` to
#' a \code{\link[grid]{unit}} object using \code{\link[ggplot2]{theme}},
#' for example: `prism.ticks.length = unit(2, "pt")`. The major tick
#' lengths are adjusted using the standard `axis.ticks.length`.
#'
#' @inheritParams ggplot2::guide_axis
#'
#' @return Returns a \emph{prism_offset_minor} guide class object.
#'
#' @example inst/examples/ex-guide_prism_offset_minor.R
#'
#' @export
guide_prism_offset_minor <- function(title = waiver(), check.overlap = FALSE,
                                     angle = NULL, n.dodge = 1, order = 0,
                                     position = waiver()) {
  if (packageVersion("ggplot2") < "3.3.0") {
    stop("ggplot2 >= 3.3.0 needed for this function.", call. = FALSE)
  }

  structure(
    list(
      title = title,

      # customizations
      check.overlap = check.overlap,
      angle = angle,
      n.dodge = n.dodge,

      # general
      order = order,
      position = position,

      # parameter
      available_aes = c("x", "y"),

      name = "axis"
    ),
    class = c("guide", "prism_offset_minor", "prism_axis", "axis")
  )
}

#' @keywords internal
#' @export
guide_train.prism_offset_minor <- function(guide, scale, aesthetic = NULL) {

  aesthetic <- aesthetic %||% scale$aesthetics[1]

  # define major and minor breaks
  major_breaks <- scale$get_breaks()
  major_breaks <- major_breaks[!is.na(major_breaks)]

  if (is.null(scale$minor_breaks)) {
    stop("No minor breaks exist, guide_prism_offset_minor needs minor breaks to work")
  }
  minor_breaks <- setdiff(scale$get_breaks_minor(), major_breaks)

  # define all breaks
  breaks <- union(major_breaks, minor_breaks)

  # indicate which breaks are major
  is_major <- breaks %in% major_breaks

  empty_ticks <- base::data.frame(
    list(aesthetic = numeric(0), .value = numeric(0), .label = character(0))
  )
  names(empty_ticks) <- c(aesthetic, ".value", ".label")

  if (length(intersect(scale$aesthetics, guide$available_aes)) == 0) {
    warn(glue(
      "axis guide needs appropriate scales: ",
      glue_collapse(guide$available_aes, ", ", last = " or ")
    ))
    guide$key <- empty_ticks
  } else if (length(breaks) == 0) {
    guide$key <- empty_ticks
  } else {
    mapped_breaks <- if (scale$is_discrete()) scale$map(breaks) else breaks

    ticks <- base::data.frame(setNames(list(mapped_breaks), aesthetic))
    ticks$.value <- breaks
    # get major break labels and make minor breaks blank
    ticks$.label <- c(scale$get_labels(major_breaks),
                      rep("", times = length(minor_breaks)))

    guide$key <- ticks[is.finite(ticks[[aesthetic]]), ]
  }

  guide$name <- paste0(guide$name, "_", aesthetic)
  guide$is_major <- is_major
  guide$hash <- digest(list(guide$title, guide$key$.value,
                            guide$key$.label, guide$name, guide$is_major))
  guide
}

#' @keywords internal
#' @export
guide_gengrob.prism_offset_minor <- function(guide, theme) {
  aesthetic <- names(guide$key)[!grepl("^\\.", names(guide$key))][1]

  draw_prism_offset_minor(
    break_positions = guide$key[[aesthetic]],
    break_labels = guide$key$.label[guide$is_major],
    breaks_major = guide$is_major,
    axis_position = guide$position,
    theme = theme,
    check.overlap = guide$check.overlap,
    angle = guide$angle,
    n.dodge = guide$n.dodge
  )
}

#' Grob for offset axes with minor ticks
#'
#' @param break_positions position of ticks
#' @param break_labels labels at ticks
#' @param breaks_major logical vector indicating major ticks versus minor ticks
#' @param axis_position position of axis (top, bottom, left or right)
#' @param theme A complete \code{\link[ggplot2]{theme}} object
#' @param check.overlap silently remove overlapping labels,
#'   (recursively) prioritizing the first, last, and middle labels.
#' @param angle Compared to setting the angle in
#'   \code{\link[ggplot2]{theme}} / `element_text`,
#'   this also uses some heuristics to automatically pick the `hjust` and
#'   `vjust` that you probably want.
#' @param n.dodge The number of rows (for vertical axes) or columns (for
#'   horizontal axes) that should be used to render the labels. This is
#'   useful for displaying labels that would otherwise overlap.
#' @keywords internal
draw_prism_offset_minor <- function(break_positions, break_labels, breaks_major,
                                    axis_position, theme,
                                    check.overlap = FALSE, angle = NULL, n.dodge = 1) {

  axis_position <- match.arg(axis_position, c("top", "bottom", "right", "left"))
  aesthetic <- if (axis_position %in% c("top", "bottom")) "x" else "y"

  # resolve elements
  line_element_name <- paste0("axis.line.", aesthetic, ".", axis_position)
  tick_element_name <- paste0("axis.ticks.", aesthetic, ".", axis_position)
  tick_length_element_name <- paste0("axis.ticks.length.", aesthetic, ".", axis_position)
  label_element_name <- paste0("axis.text.", aesthetic, ".", axis_position)
  prism_tick_element_name <- paste0("prism.ticks.", aesthetic, ".", axis_position)
  prism_tick_length_element_name <- paste0("prism.ticks.length.", aesthetic, ".", axis_position)

  line_element <- calc_element(line_element_name, theme)
  tick_element <- calc_element(tick_element_name, theme)
  tick_length <- calc_element(tick_length_element_name, theme)
  label_element <- calc_element(label_element_name, theme)
  prism_tick_element <- calc_element(prism_tick_element_name, theme)
  prism_tick_length <- calc_element(prism_tick_length_element_name, theme)

  # override label element parameters for rotation
  if (inherits(label_element, "element_text")) {
    label_overrides <- axis_label_element_overrides(axis_position, angle)
    # label_overrides is an element_text, but label_element may not be;
    # to merge the two elements, we just copy angle, hjust, and vjust
    # unless their values are NULL
    if (!is.null(label_overrides$angle)) {
      label_element$angle <- label_overrides$angle
    }
    if (!is.null(label_overrides$hjust)) {
      label_element$hjust <- label_overrides$hjust
    }
    if (!is.null(label_overrides$vjust)) {
      label_element$vjust <- label_overrides$vjust
    }
  }

  # conditionally set parameters that depend on axis orientation
  is_vertical <- axis_position %in% c("left",  "right")

  position_dim <- if (is_vertical) "y" else "x"
  non_position_dim <- if (is_vertical) "x" else "y"
  position_size <- if (is_vertical) "height" else "width"
  non_position_size <- if (is_vertical) "width" else "height"
  gtable_element <- if (is_vertical) gtable_row else gtable_col
  measure_gtable <- if (is_vertical) gtable_width else gtable_height
  measure_labels_non_pos <- if (is_vertical) grobWidth else grobHeight

  # conditionally set parameters that depend on which side of the panel
  # the axis is on
  is_second <- axis_position %in% c("right", "top")

  tick_direction <- if (is_second) 1 else -1
  non_position_panel <- if (is_second) unit(0, "npc") else unit(1, "npc")
  tick_coordinate_order <- if (is_second) c(2, 1) else c(1, 2)

  # conditionally set the gtable ordering
  labels_first_gtable <- axis_position %in% c("left", "top") # refers to position in gtable

  # set common parameters
  n_breaks <- length(break_positions[breaks_major])
  n_minor_breaks <- length(break_positions[!breaks_major])
  opposite_positions <- c("top" = "bottom", "bottom" = "top", "right" = "left", "left" = "right")
  axis_position_opposite <- unname(opposite_positions[axis_position])

  # draw elements
  line_grob <- exec(
    element_grob, line_element,
    !!position_dim := unit(c(min(break_positions),
                             max(break_positions)), "npc"),
    !!non_position_dim := unit.c(non_position_panel, non_position_panel)
  )

  if (n_breaks == 0) {
    return(
      absoluteGrob(
        gList(line_grob),
        width = grobWidth(line_grob),
        height = grobHeight(line_grob)
      )
    )
  }

  # break_labels can be a list() of language objects
  if (is.list(break_labels)) {
    if (any(vapply(break_labels, is.language, logical(1)))) {
      break_labels <- do.call(expression, break_labels)
    } else {
      break_labels <- unlist(break_labels)
    }
  }

  # calculate multiple rows/columns of labels (which is usually 1)
  dodge_pos <- rep(seq_len(n.dodge), length.out = n_breaks)
  dodge_indices <- split(seq_len(n_breaks), dodge_pos)

  label_grobs <- lapply(dodge_indices, function(indices) {
    draw_axis_labels(
      break_positions = break_positions[breaks_major][indices],
      break_labels = break_labels[indices],
      label_element = label_element,
      is_vertical = is_vertical,
      check.overlap = check.overlap
    )
  })

  ticks_grob <- exec(
    element_grob, tick_element,
    !!position_dim := rep(unit(break_positions, "native"), each = 2),
    !!non_position_dim := unit.c(
      rep(
        unit.c(non_position_panel + (tick_direction * tick_length), non_position_panel)[tick_coordinate_order],
        times = n_breaks
      ),
      rep(
        # let minor ticks have a different length from major ticks
        unit.c(non_position_panel + (tick_direction * prism_tick_length), non_position_panel)[tick_coordinate_order],
        times = n_minor_breaks
      )
    ),
    id.lengths = rep(2, times = length(break_positions))
  )

  # create gtable
  non_position_sizes <- paste0(non_position_size, "s")
  label_dims <- do.call(unit.c, lapply(label_grobs, measure_labels_non_pos))
  grobs <- c(list(ticks_grob), label_grobs)
  grob_dims <- unit.c(tick_length, label_dims)

  if (labels_first_gtable) {
    grobs <- rev(grobs)
    grob_dims <- rev(grob_dims)
  }

  gt <- exec(
    gtable_element,
    name = "axis",
    grobs = grobs,
    !!non_position_sizes := grob_dims,
    !!position_size := unit(1, "npc")
  )

  # create viewport
  justvp <- exec(
    viewport,
    !!non_position_dim := non_position_panel,
    !!non_position_size := measure_gtable(gt),
    just = axis_position_opposite
  )

  absoluteGrob(
    gList(line_grob, gt),
    width = gtable_width(gt),
    height = gtable_height(gt),
    vp = justvp
  )
}

#' Prism themes
#'
#' A collection of ggplot2 themes that use palettes which mirror the
#' colour schemes available in GraphPad Prism.
#'
#' @param palette `string`. Palette name, use
#' `names(ggprism_data$themes)` to show all valid palette names.
#' @param base_size `numeric`. Base font size, given in `"pt"`.
#' @param base_family `string`. Base font family, default is `"sans"`.
#' @param base_fontface `string`. Base font face, default is `"bold"`.
#' @param base_line_size `numeric`. Base linewidth for line elements
#' @param base_rect_size `numeric`. Base linewidth for rect elements
#' @param axis_text_angle `integer`. Angle of axis text in degrees.
#' One of: `0, 45, 90, 270`.
#' @param border `logical`. Should a border be drawn around the plot?
#' Clipping will occur unless e.g. `coord_cartesian(clip = "off")` is used.
#'
#' @return Returns a list-like object of class _theme_.
#'
#' @noRd
#'
#' @example inst/examples/ex-theme_prism.R
#'
#' @export
theme_prism_scCustomize <- function(palette = "black_and_white", base_size = 14,
                                    base_family = "sans", base_fontface = "bold",
                                    base_line_size = base_size / 14,
                                    base_rect_size = base_size / 14,
                                    axis_text_angle = 0,
                                    border = FALSE) {

  # Ensure x axis text is at a sensible angle
  angle <- axis_text_angle[1]
  if(!angle %in% c(0, 45, 90, 270))
    stop(sprintf("'axis_text_angle' must be one of [%s]",
                 paste(c(0, 45, 90, 270), collapse=", ")),
         ".\nFor other angles, use the guide_axis() function in ggplot2 instead",
         call. = FALSE)

  # Get element colours from palette
  if (!palette %in% names(ggprism::ggprism_data$themes)) {
    stop("The palette ", paste(palette), " does not exist.
         See names(ggprism_data$themes) for valid palette names")
  }
  colours <- deframe(ggprism::ggprism_data$themes[[palette]])

  # Draw border or not
  if(!is_bool(border)) {
    stop("border must be either: TRUE or FALSE")
  } else {
    if(border){
      panel.border <- element_rect(fill = NA)
      axis.line <- element_blank()
    }
    else if (!border) {
      panel.border <- element_blank()
      axis.line <- element_line()
    }
  }

  t <- theme(
    # Base elements (to be inherited by other elements)
    line = element_line(colour = colours["axisColor"], size = base_line_size,
                        linetype = 1, lineend = "square"),
    rect = element_rect(fill = "white", colour = colours["axisColor"],
                        size = base_rect_size, linetype = 1),
    text = element_text(family = base_family, face = base_fontface,
                        colour = colours["graphTitleColor"], size = base_size,
                        lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0,
                        margin = margin(), debug = FALSE),

    # Prism custom theme elements
    prism.ticks.length = unit(base_size / 5, "pt"),

    # Normal ggplot2 theme elements
    axis.line =          axis.line,
    axis.line.x =        NULL,
    axis.line.y =        NULL,
    axis.text =          element_text(size = rel(0.95), colour = colours["axisLabelColor"]),
    axis.text.x =        element_text(margin = margin(t = 0.8 * base_size / 4),
                                      angle = axis_text_angle,
                                      hjust = ifelse(axis_text_angle %in% c(45, 90, 270), 1, 0.5),
                                      vjust = ifelse(axis_text_angle %in% c(0, 90, 270), 0.5, 1)),
    axis.text.x.top =    element_text(margin = margin(b = 0.8 * base_size / 4), vjust = 0),
    axis.text.y =        element_text(margin = margin(r = 0.5 * base_size / 4), hjust = 1),
    axis.text.y.right =  element_text(margin = margin(l = 0.5 * base_size / 4), hjust = 0),
    axis.ticks =         element_line(),
    axis.ticks.length =  unit(base_size / 2.5, "pt"),
    axis.ticks.length.x = NULL,
    axis.ticks.length.x.top = NULL,
    axis.ticks.length.x.bottom = NULL,
    axis.ticks.length.y = NULL,
    axis.ticks.length.y.left = NULL,
    axis.ticks.length.y.right = NULL,
    axis.title =         element_text(colour = colours["axisTitleColor"]),
    axis.title.x =       element_text(margin = margin(t = base_size * 0.6),
                                      vjust = 1),
    axis.title.x.top =   element_text(margin = margin(b = base_size * 0.6),
                                      vjust = 0),
    axis.title.y =       element_text(angle = 90,
                                      margin = margin(r = base_size * 0.6),
                                      vjust = 1),
    axis.title.y.right = element_text(angle = -90,
                                      margin = margin(l = base_size * 0.6),
                                      vjust = 0),

    legend.background =  element_blank(),
    legend.spacing =     unit(base_size, "pt"),
    legend.spacing.x =   NULL,
    legend.spacing.y =   NULL,
    legend.margin =      margin(base_size / 2,  base_size / 2,
                                base_size / 2, base_size / 2),
    legend.key =         element_blank(),
    legend.key.size =    unit(1.2, "lines"),
    legend.key.height =  NULL,
    legend.key.width =   unit(base_size * 1.8, "pt"),
    legend.text =        element_text(size = rel(0.8), face = "plain"),
    legend.text.align =  NULL,
    legend.title =       element_blank(),
    legend.title.align = NULL,
    legend.position =    "right",
    legend.direction =   NULL,
    legend.justification = "center",
    legend.box =         NULL,
    legend.box.margin =  margin(0, 0, 0, 0, "cm"),
    legend.box.background = element_blank(),
    legend.box.spacing = unit(base_size, "pt"),

    panel.background = element_rect(fill = ifelse(palette == "office", colours["plottingAreaColor"], NA),
                                    colour = NA),
    panel.border =       panel.border,
    panel.grid =         element_blank(),
    panel.grid.minor =   element_blank(),
    panel.spacing =      unit(base_size / 2, "pt"),
    panel.spacing.x =    NULL,
    panel.spacing.y =    NULL,
    panel.ontop    =     FALSE,

    strip.background =   element_blank(),
    strip.text =         element_text(colour = colours["axisTitleColor"],
                                      size = rel(0.8),
                                      margin = margin(base_size / 2.5, base_size / 2.5,
                                                      base_size / 2.5, base_size / 2.5)),
    strip.text.x =       element_text(margin = margin(b = base_size / 3)),
    strip.text.y =       element_text(angle = -90, margin = margin(l = base_size / 3)),
    strip.text.y.left =  element_text(angle = 90),
    strip.placement =    "inside",
    strip.placement.x =  NULL,
    strip.placement.y =  NULL,
    strip.switch.pad.grid = unit(base_size / 4, "pt"),
    strip.switch.pad.wrap = unit(base_size / 4, "pt"),

    plot.background =    element_rect(fill = colours["pageBackgroundColor"],
                                      colour = NA),
    plot.title =         element_text(size = rel(1.2),
                                      hjust = 0.5, vjust = 1,
                                      margin = margin(b = base_size)),
    plot.title.position = "panel",
    plot.subtitle =      element_text(hjust = 0.5, vjust = 1,
                                      margin = margin(b = base_size / 2)),
    plot.caption =       element_text(size = rel(0.8),
                                      hjust = 1, vjust = 1,
                                      margin = margin(t = base_size / 2)),
    plot.caption.position = "panel",
    plot.tag =           element_text(size = rel(1.2),
                                      hjust = 0.5, vjust = 0.5),
    plot.tag.position =  'topleft',
    plot.margin =        margin(base_size / 2, base_size / 2,
                                base_size / 2, base_size / 2),

    complete = TRUE
  )

  parent <- ggprism::ggprism_data$themes[["all_null"]]
  if (!"legend.text.align" %in% rlang::fn_fmls_names(theme)) {
    t$legend.text.align  <- parent$legend.text.align <- NULL
    t$legend.title.align <- parent$legend.title.align <- NULL
  }

  # make sure all elements are set to NULL if not explicitly defined
  parent %+replace% t
}
