## <<<<<<<<<<<<<<< HEAD

# load packages
library(tidyverse)
library(ggpubr)
library(readxl)
library(ggh4x)
library(ggsci)

dat <- read.csv("Data/reference_solar_spectrum.csv")

x <- dat$lambda
y <- dat$terrestrial_direct

UV <- integrate(approxfun(x,y), 280, 400, subdivisions = 220)$value
VIS <-integrate(approxfun(x,y), 400, 780, subdivisions = 400)$value
IR <- integrate(approxfun(x,y), 780, 4000, subdivisions = 1500)$value

total <- UV+VIS+IR

UV/total*100
VIS/total*100
IR/total*100




lambda_to_eV <- function(lambda_nm){
  lambda <- lambda_nm * 1e-9
  h <- 6.6261e-34 # J*s
  c <- 299792458 # m/s
  E <- (h*c)/lambda
  eV <- 1.602176565e-19 # J
  E/eV
}

svg("Results/spectral_irradiance.svg",  width=7, height=6)
dat %>% ggplot() + 
  geom_line(data = dat, aes(lambda, terrestrial), color = "grey") +
  geom_line(data = dat, aes(lambda, extraterrestrial), color = "grey") +
  geom_ribbon(aes(x=lambda, ymax=extraterrestrial, ymin=terrestrial), fill = "yellow", alpha=.5) +
  geom_ribbon(aes(x=lambda, ymax=terrestrial, ymin=0), fill = "red",  alpha=.5) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0,3000,250), limits = c(250,3100), expand = c(0,0),  guide = "axis_minor",
                     sec.axis=sec_axis(~ lambda_to_eV(.), name="Energie (eV)", 
                                       breaks = round(lambda_to_eV(c(seq(250,3000,250))), 2))) +
  scale_y_continuous(breaks = seq(0,2, 0.2), expand = c(0,0),  guide = "axis_minor",
                     sec.axis=dup_axis(name = "")) +
  xlab(paste("Wellenlänge", "\u03bb", "(nm)")) +
  ylab(expression(Spektrale ~ Bestrahlungsstärke ~ (W ~ m^-2 ~ nm^-1))) +
  geom_vline(xintercept = 400, lty = 2, size = 0.8) +
  geom_vline(xintercept = 780, lty = 2, size = 0.8) +
  theme(text = element_text(size = 14)) +
  annotate(geom = "text", x = 330, y = 2.1, label = "UV", size = 5) +
  annotate(geom = "text", x = 600, y = 2.1, label = "VIS", size = 5) +
  annotate(geom = "text", x = 950, y = 2.1, label = "Infrarot", size = 5)
dev.off()



 
svg("Results/spectral_irradiance_legend.svg",  width=7, height=6)
dat %>% 
  gather("type", "irradiance", -lambda) %>%
  mutate(type = ifelse(type == "terrestrial", "Sonnenstrahlung auf der Erdoberfläche", "Extraterrestrische Sonnenstrahlung")) %>%
  ggplot(aes(lambda, irradiance, fill = type)) +
  geom_bar(stat = "identity", alpha = 0.5) + scale_fill_manual(values = c("yellow", "red")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position = c(0.8, 0.8)) 
dev.off()
  

