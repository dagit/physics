module Numeric.Physics.Particle where

import Numeric.LinearAlgebra.Vector

data Particle a = Particle
  { pPos     :: !(Vec3 a) -- ^ Position of particle in the world
  , pVel     :: !(Vec3 a) -- ^ Velocity of particle
  , pAcc     :: !(Vec3 a) -- ^ Acceleration of particle
  , pDamp    :: !a        -- ^ damping factor for this particle
  , pInvMass :: !a        -- ^ inverse mass (1/m) of this particle
  } deriving (Read,Show,Ord,Eq)

-- | Takes a particle, a force, and a time duration
-- and gives the particle updated for time slice
integrate :: Floating a => Particle a -> Vec3 a -> a -> Particle a
integrate p f dt = 
  let pos  = (dt         *> pVel p) <+> pPos p
      acc' = (pInvMass p *> f)      <+> pAcc p
      vel  = (dt         *> acc')   <+> pVel p
  in p { pPos = pos
       , pVel = ((pDamp p)**dt) *> vel
       }
{-# SPECIALIZE INLINE integrate :: Particle Double -> Vec3 Double -> Double -> Particle Double #-}
{-# SPECIALIZE INLINE integrate :: Particle Float  -> Vec3 Float  -> Float  -> Particle Float  #-}
